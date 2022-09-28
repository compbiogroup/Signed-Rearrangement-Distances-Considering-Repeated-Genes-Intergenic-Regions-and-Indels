{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

-- |
-- Module      : MCISP2K
-- Description : Heuristic for Genome partition problems.
-- Copyright   : (c) Gabriel Siqueira, 2021
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
module MCISP2K where

import Balancer (reduceGenesForDeletion)
import Control.Concurrent.ParallelIO.Global (parallel_, stopGlobalPool)
import Control.DeepSeq (force)
import qualified Data.ByteString.Char8 as BS
import Data.Time (diffUTCTime, getCurrentTime)
import Genomes (Genome, Sign(..), readGenome, writeGenome, zeroIR)
import LocalBase
import Options.Applicative
import Partition (getPartition, reduced, mapToPerm, PartitionType(..))
import Text.Printf (printf)

data Args = Args
  { input :: String,
    output :: String,
    noParallel :: Bool,
    signed :: Sign,
    partType :: PartitionType,
    asPerm :: Bool
  }

argsParser :: Parser Args
argsParser =
  Args
    <$> strOption
      ( long "input"
          <> short 'i'
          <> metavar "IFILE"
          <> help "Input file. Each 4 lines of the input file correspond to a instance, each line has a list of comma or space separated values, and represent in order the origin string, the origin intergenic region list, the target string, and the target intergenic region list."
      )
      <*> strOption
        ( long "outfile"
            <> short 'o'
            <> metavar "OFILE"
            <> help "Output file. For each instance five lines are produces in the file, the first four lines correspond to two reduced genomes produced with the partition algorithm (each character of the string correspond to a block of the partition and each integer of the intergenic regions list correspond to a breakpoint). The last line shows the wall clock time required to produce the partition."
        )
      <*> switch
        ( long "no-par"
            <> help "Do not process the genomes in parallel."
        )
      <*> flag Unsigned Signed 
        ( long "signed"
            <> short 's'
            <> help "Whether the input Strings are signed."
        )
      <*> flag MCISP RMCISP
        ( long "rev"
            <> help "Whether to use reverse partition."
        )
      <*> switch
        ( long "perm"
            <> help "Whether to produce permutations from a mapping of the original strings instead of reduced genomes."
        )

opts :: ParserInfo Args
opts =
  info
    (argsParser <**> helper)
    ( fullDesc
        <> progDesc
          "Algorithm for genome partition problems."
    )

main :: IO ()
main = do
  args <- execParser opts
  contents <- BS.readFile (input args)
  let pairs = zip [(1::Int)..] . toQuadruples . filter ((/= '#') . BS.head) . BS.lines $ contents
  if noParallel args
    then mapM_ (runOne args) pairs
    else do
        parallel_ $ map (runOne args) pairs
        stopGlobalPool
  where
    runOne args (i, bstrs) = do
      start <- getCurrentTime
      let !bstrs' = force $ simplifyGenomes (asPerm args) (partType args) (signed args) bstrs
      end <- getCurrentTime
      let time = BS.pack . show . realToFrac $ diffUTCTime end start
      BS.writeFile (output args ++ "_" ++ printf "%04d" i) . BS.unlines $ fromAns (bstrs', "# Time: " <> (BS.pack . show $ time))

    toQuadruples (s1 : i1 : s2 : i2 : ss) = (s1, i1, s2, i2) : toQuadruples ss
    toQuadruples [] = []
    toQuadruples _ = error "Incorrect number of lines."

    fromAns ((s1, i1, s2, i2), time)  = s1 : i1 : s2 : i2 : [time]

simplifyGenomes :: Bool -> PartitionType -> Sign -> (BS.ByteString, BS.ByteString, BS.ByteString, BS.ByteString) -> (BS.ByteString, BS.ByteString, BS.ByteString, BS.ByteString)
simplifyGenomes asPerm ptype signed (s1, i1, s2, i2) = (s1', i1', s2', i2')
  where
    (s1', i1') = writeGenome False g'
    (s2', i2') = writeGenome False h'
    (g_bal, h_bal) = reduceGenesForDeletion g h
    -- g_bal_z = zeroIR g_bal
    -- h_bal_z = zeroIR h_bal
    -- part = getPartition ptype g_bal_z h_bal_z
    part = getPartition ptype g_bal h_bal
    (g', h') = if asPerm then mapToPerm g_bal h_bal part else reduced part
    g = readGenome True signed s1 i1
    h = readGenome True signed s2 i2
