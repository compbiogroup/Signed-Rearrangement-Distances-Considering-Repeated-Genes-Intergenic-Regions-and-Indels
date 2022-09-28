{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TupleSections #-}

module DB where

import Control.Monad (replicateM, (>=>))
import Control.Monad.Random (Rand, StdGen, evalRandIO, getRandomR, getRandoms)
import Data.ByteString (ByteString)
import qualified Data.ByteString.Char8 as BS
import Data.Coerce (coerce)
import Data.List (sort, unfoldr)
import Genomes
import LocalBase
import Options.Applicative

opts :: ParserInfo Args
opts =
  info
    (argParser <**> helper)
    ( fullDesc
        <> progDesc "Generate database with genomes. Used for tests of solutions for rearrangement problems."
    )

data Parameters = DB1 RepDB | DB2 RandDB

data RepDB = RepDB
  { db_low :: Int,
    db_high :: Int,
    db_rep :: Int
  }

newtype RandDB = RandDB {db_lim :: Int}

data Args = Args
  { db_par :: Parameters,
    db_num_pairs :: Int,
    db_size :: Size,
    db_nop :: Int,
    db_porc :: Int,
    db_sign :: Sign,
    db_indel :: Int,
    db_zeros :: Bool,
    db_output :: String
  }

repDBParser :: Parser RepDB
repDBParser =
  RepDB
    <$> option
      auto
      ( long "low"
          <> short 'l'
          <> metavar "L"
          <> help "Minimum number of replicas."
      )
    <*> option
      auto
      ( long "high"
          <> short 'h'
          <> metavar "H"
          <> help "Maximum number of replicas."
      )
    <*> option
      auto
      ( long "replicas"
          <> short 'r'
          <> metavar "REP"
          <> help "Number of genes having more than one occurrence."
      )

randDBParser :: Parser RandDB
randDBParser =
  RandDB
    <$> option
      auto
      ( long "alphabet-size"
          <> short 'l'
          <> metavar "ALPH"
          <> help "Alphabet size."
      )

argParser :: Parser Args
argParser =
  Args
    <$> subparser
      ( command
          "RepDB"
          ( info
              (DB1 <$> repDBParser <**> helper)
              (progDesc "DB with fix number of replicas.")
          )
          <> command
            "RandDB"
            ( info
                (DB2 <$> randDBParser <**> helper)
                (progDesc "DB with random generation.")
            )
      )
      <*> option
        auto
        ( long "number_genomes"
            <> short 'k'
            <> metavar "K"
            <> help "Number genome pairs to generate."
        )
      <*> option
        auto
        ( long "size_genome"
            <> short 'n'
            <> metavar "N"
            <> help "Size of the genomes."
        )
      <*> option
        auto
        ( long "number_op"
            <> short 'r'
            <> metavar "R"
            <> help "Number of operations to apply (-1 to use a random list)."
        )
      <*> option
        auto
        ( long "porcentage_rev"
            <> short 'p'
            <> metavar "P"
            <> showDefault
            <> value 100
            <> help "Porcentage of reversions in the operations."
        )
      <*> flag
        Unsigned
        Signed
        ( long "signed"
            <> short 's'
            <> help "Whether the input Strings are signed."
        )
      <*> option
        auto
        ( long "indel"
            <> short 'd'
            <> metavar "D"
            <> help "How much indels to apply (D deletions follow by D insertions)."
        )
      <*> switch
        ( long "zeros"
            <> short 'z'
            <> help "Whether to produce intergenic regions with zeros."
        )
      <*> strOption
        ( long "outfile"
            <> short 'o'
            <> metavar "oFILE"
            <> help "Output file"
        )

main :: IO ()
main = do
  args <- execParser opts
  db <- evalRandIO . fmap (BS.unlines . fromQuadruples) $ replicateM (db_num_pairs args) (genPair args)
  BS.writeFile (db_output args) db
  where
    fromQuadruples ((s1, i1, s2, i2) : ss) = s1 : i1 : s2 : i2 : fromQuadruples ss
    fromQuadruples [] = []

genPair :: Args -> Rand StdGen (ByteString, ByteString, ByteString, ByteString)
genPair Args {..} = do
  g <- case db_par of
    (DB1 (RepDB l h d)) -> randomGenomeWithReplicas db_zeros db_size d l h db_sign
    (DB2 (RandDB lim)) -> randomGenome db_zeros db_size lim db_sign
  h <-
    applyIndels
      =<< if db_nop == -1
        then rearrangeGenome g
        else applyOperations g
  let (s1, i1) = writeGenome True g
  let (s2, i2) = writeGenome True h
  return (s1, i1, s2, i2)
  where
    r_r = (db_nop * db_porc) `div` 100
    r_t = db_nop - r_r

    applyIndels g = do
      let dels = unfoldr dels_for_one db_indel
      g' <- foldr (=<<) (return g) dels
      let ins = unfoldr ins_for_one (db_indel, (+1) . (!! 1) . reverse . sort . alphabet $ g)
      foldr (=<<) (return g') ins
    dels_for_one 0 = Nothing
    dels_for_one d = Just . (,d -1) $ \g -> do
      i <- getRandomR (2 :: Idx, coerce $ genomeSize g - 1)
      ir <- getRandomR (0, irByIdx g (i -1) + irByIdx g i)
      return $ deletion i (i + 1) ir g
    ins_for_one (0,_) = Nothing
    ins_for_one (d,next) = Just . (,(d-1,next+1)) $ \g -> do
      i <- getRandomR (1 :: Idx, coerce $ genomeSize g - 1)
      ir1 <- if db_zeros then return 0 else getRandomR (0, 100)
      ir2 <- if db_zeros then return 0 else getRandomR (max 0 (irByIdx g i - ir1), 100)
      return $ insertion i (fromLists False (genomeIsSigned g) [coerce next] [ir1,ir2]) g

    applyOperations :: Genome -> Rand StdGen Genome
    applyOperations g = do
      coins <- getRandoms
      let ops = unfoldr operations_for_one (r_t, r_r, coins)
      foldr (=<<) (return g) ops
    operations_for_one :: (Int, Int, [Bool]) -> Maybe (Genome -> Rand StdGen Genome, (Int, Int, [Bool]))
    operations_for_one (_, _, []) = Nothing
    operations_for_one (r_t', r_r', coin : coins)
      | r_t' == 0 && r_r' == 0 = Nothing
      | r_t' == 0 || r_r' /= 0 && coin = Just . (,(r_t', r_r' - 1, coins)) $ \g -> do
        i <- getRandomR (2 :: Idx, coerce db_size - 2)
        j <- getRandomR (i + 1, coerce db_size - 1)
        x <- getRandomR (0, irByIdx g (i -1))
        y <- getRandomR (0, irByIdx g j)
        return $ reversal i j x y g
      | otherwise = Just . (,(r_t' - 1, r_r', coins)) $ \g -> do
        i <- getRandomR (2 :: Idx, coerce db_size - 2)
        j <- getRandomR (i + 1, coerce db_size - 1)
        k <- getRandomR (j + 1, coerce db_size)
        x <- getRandomR (0, irByIdx g (i -1))
        y <- getRandomR (0, irByIdx g (j -1))
        z <- getRandomR (0, irByIdx g (k -1))
        return $ transposition i j k x y z g
