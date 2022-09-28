{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Module      : Genomes
-- Description : Representation of a genome
-- Copyright   : (c) Gabriel Siqueira, 2021
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
--
-- A genome comprises of a gene list and a integer list (correspondent to the sizes of intergenic regions).
module Genomes
  ( Duo,
    Gene,
    Genome (genomeIsSigned),
    IR (..),
    Idx (..),
    Sign (..),
    allSubGenomes,
    balanced,
    breakGenome,
    combineGenomes,
    combineGenomesL,
    duoByIdx,
    irByIdx,
    duosList,
    duoIdx,
    duoIr,
    duoToBS,
    estimateITD,
    zeroIR,
    fromLists,
    genePair,
    genomeSize,
    alphabet,
    occurenceMap,
    occurenceMapLookup,
    occurenceMapAdjust,
    intToGene,
    interleaveListRepresentation,
    interleaveListToGenome,
    maybeEstimateITD,
    occurenceMax,
    geneMaxValue,
    randomGenome,
    randomGenomeWithReplicas,
    rearrangeGenome,
    readGenome,
    writeGenome,
    sliceGenome,
    reversal,
    transposition,
    insertion,
    deletion,
    subGenCount,
    subGenome,
    toLists,
    validBeginIR,
    validEndIR,
    -- weigth,
    subGenomeFind,
    GenomeMap,
    GeneListForMap(..),
    gmEmpty,
    gmEmptyRev,
    gmLookup,
    gmAlter,
    gmInsert,
    gmLookupInsert,
  )
where

import Control.Exception (assert)
import Control.Monad.Random (MonadRandom, getRandomRs, getRandoms)
import Data.ByteString.Builder (intDec, toLazyByteString)
import qualified Data.ByteString.Char8 as BS
import qualified Data.ByteString.Lazy as LBS
import Data.Char (isNumber)
import Data.Coerce (coerce)
import Data.Foldable (foldl', toList)
import Data.Hashable (Hashable)
import Data.IntMap.Strict (IntMap)
import qualified Data.IntMap.Strict as IntMap
import qualified Data.List as List
import Data.Maybe (fromJust)
import Data.Text (Text)
import qualified Data.Text as Text
import Data.Vector (Vector, (!))
import qualified Data.Vector as Vec
import qualified Data.Vector.Mutable as MVec
import GHC.Generics (Generic)
import LocalBase
import System.Random (Random)
import System.Random.Shuffle (shuffleM)

newtype IR = IR Int deriving newtype (Eq, Show, Read, Num, Ord, Random)

newtype Idx = Idx Int deriving newtype (Eq, Show, Read, Num, Ord, Enum, Random, Integral, Real)

newtype Gene = Gene Int deriving newtype (Eq, Show, Read, Hashable, Ord, Num, Random, Bounded, Enum)

type Gstring = Vector Gene

type IRList = Vector IR

data Sign = Signed | Unsigned deriving (Eq, Show, Enum)

-- | Representation of a genome the gstring must be a non empty sequence of genes and
--  size of irList must be size of gstring minus one
data Genome = Genome
  { gstring :: Gstring,
    irList :: IRList,
    interleaveListRepresentation :: Vector BS.ByteString,
    genomeIsSigned :: Sign
  }
  deriving (Eq)

data Duo = Duo
  { genePair :: (Gene, Gene),
    duoIr :: IR,
    duoIdx :: Idx,
    parGenomeSize :: Size,
    duoIsSigned :: Sign
  }
  deriving (Show)

----------------------------------------
--     Instantiations of Type Classes --
----------------------------------------

instance Orientable Gene where
  getOri a = if a >= 0 then LR else RL
  invOri a = - a

instance Orientable Duo where
  getOri duo =
    case duoIsSigned duo of
      Signed
        | getOri a1 == LR && getOri a2 == LR -> LR
        | getOri a1 == RL && getOri a2 == RL -> RL
        | otherwise -> if canonicOri a1 < canonicOri a2 then LR else RL
      Unsigned -> if canonicOri a1 < canonicOri a2 then LR else RL
    where
      (a1, a2) = genePair duo

  invOri duo@Duo {..} = case duoIsSigned of
    Signed -> duo {genePair = (invOri a2, invOri a1), duoIdx = idx'}
    Unsigned -> duo {genePair = (a2, a1), duoIdx = idx'}
    where
      (a1, a2) = genePair
      idx' = coerce parGenomeSize - duoIdx

instance Semigroup Sign where
  Signed <> Signed = Signed
  Signed <> Unsigned = Signed
  Unsigned <> Signed = Signed
  Unsigned <> Unsigned = Unsigned

instance Show Genome where
  show g =
    unwords . (("(" ++ head str_s ++ ")") :) $
      zipWith (\ir a -> "- " ++ ir ++ " - (" ++ a ++ ")") str_i (tail str_s)
    where
      str_s = Vec.toList $ (\i -> if i == maxBound then "inf" else show i) <$> gstring g
      str_i = Vec.toList $ show <$> irList g

instance Orientable Genome where
  getOri g =
    if
        | gstring g < gstring rg ->
          case genomeIsSigned g of
            Signed -> RL
            Unsigned -> LR
        | gstring g > gstring rg ->
          case genomeIsSigned g of
            Signed -> LR
            Unsigned -> RL
        | otherwise -> if irList g <= irList rg then LR else RL
    where
      rg = invOri g
  invOri g = Genome rs ri ilr (genomeIsSigned g)
    where
      ilr = Vec.fromList $ interleavelists ls_bs li_bs
      ls_bs = fmap geneToBS . toList $ rs
      li_bs = fmap irToBS . toList $ ri
      ri = Vec.reverse . irList $ g
      rs =
        ( case genomeIsSigned g of
            Signed -> fmap invOri
            Unsigned -> id
        )
          . Vec.reverse
          . gstring
          $ g

instance Orientable (Maybe Genome) where
    getOri Nothing = LR
    getOri (Just x) = getOri x
    invOri = fmap invOri

--------------------------------
--      Random Generation     --
--------------------------------

randomGenome :: MonadRandom mon => Bool -> Size -> Int -> Sign -> mon Genome
randomGenome zeros size lim signed = do
  coins <- getRandoms
  ls <- case signed of
    Unsigned -> take n <$> getRandomRs (1, coerce lim)
    Signed -> zipWith swaps coins . take n <$> getRandomRs (1, coerce lim)
  li <- take (n + 1) <$> if zeros then return (repeat 0) else getRandomRs (1, 100)
  return $ fromLists True signed ls li
  where
    n = fromIntegral size
    swaps b v = if b then v else invOri v

randomGenomeWithReplicas :: MonadRandom mon => Bool -> Size -> Int -> Int -> Int -> Sign -> mon Genome
randomGenomeWithReplicas zeros size rep low high signed = do
  coins <- getRandoms
  occs <- getRandomRs (fromIntegral low, fromIntegral high)
  let ls =
        ( case signed of
            Unsigned -> id
            Signed -> zipWith swaps coins
        )
          . (\l -> l ++ coerce [length l + 1 .. n])
          . concat
          . zipWith replicate occs
          . coerce
          $ [1 .. rep]
  li <- take (n + 1) <$> if zeros then return (repeat 0) else getRandomRs (1, 100)
  return $ fromLists True signed ls li
  where
    n = fromIntegral size
    swaps b v = if b then v else invOri v

rearrangeGenome :: MonadRandom mon => Genome -> mon Genome
rearrangeGenome g = do
  coins <- getRandoms
  let (sign, ls, li) = toLists True g
      s = sum li
  ls' <- shuffleM ls
  ls' <- case sign of
    Unsigned -> shuffleM ls
    Signed -> zipWith swaps coins <$> shuffleM ls
  x <-
    List.sort . map IR . take (coerce $ genomeSize g) <$> getRandomRs (0, coerce s)
  let li' = zipWith (-) (x ++ [s]) (0 : x)
  return $ fromLists True sign ls' li'
  where
    swaps b v = if b then v else invOri v

--------------------------
--      Conversions     --
--------------------------

readGenome :: Bool -> Sign -> BS.ByteString -> BS.ByteString -> Genome
readGenome extend sign s i =
  fromLists
    extend
    sign
    (map (Gene . readInt) . BS.splitWith (\x -> x == ',' || x == ' ') $ s)
    (map (IR . readInt) . BS.splitWith (\x -> x == ',' || x == ' ') $ i)
  where
    readInt = coerce . fst . fromJust . BS.readInt

writeGenome :: Bool -> Genome -> (BS.ByteString, BS.ByteString)
writeGenome rext g =
  ( BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec . coerce) $ ls,
    BS.unwords . fmap (LBS.toStrict . toLazyByteString . intDec . coerce) $ li
  )
  where
    (_, ls, li) = toLists rext g

duoToBS :: Duo -> [BS.ByteString]
duoToBS Duo {..} =
  let (al, ar) = genePair
   in [geneToBS al, irToBS duoIr, geneToBS ar]

intToGene :: Int -> Gene
intToGene = coerce

geneToBS :: Gene -> BS.ByteString
geneToBS = LBS.toStrict . toLazyByteString . (<>) "g" . intDec . coerce

irToBS :: IR -> BS.ByteString
irToBS = LBS.toStrict . toLazyByteString . (<>) "i" . intDec . coerce

fromLists :: Bool -> Sign -> [Gene] -> [IR] -> Genome
fromLists extend sign ls_ li = Genome (Vec.fromList ls) (Vec.fromList li) ilr sign
  where
    ilr = Vec.fromList $ interleavelists ls_bs li_bs
    ls_bs = fmap geneToBS ls
    li_bs = fmap irToBS li

    ls__ = case sign of Signed -> map canonicOri ls_; Unsigned -> ls_
    ls =
      if extend
        then 0 : (ls_ ++ [maxBound])
        else ls_

toLists :: Bool -> Genome -> (Sign, [Gene], [IR])
toLists rext g = (genomeIsSigned g, Vec.toList . (if rext then Vec.slice 1 (coerce $ genomeSize g - 2) else id) $ gstring g, Vec.toList $ irList g)

duosList :: Genome -> [Duo]
duosList g = foldr toDuo [] . zip [1 ..] . lPairs . Vec.toList $ gstring g
  where
    irs = irList g
    toDuo (i, (al, ar)) l = Duo (al, ar) (irs ! (i -1)) (Idx i) (genomeSize g) (genomeIsSigned g) : l

interleaveListToGenome :: Vector BS.ByteString -> Sign -> Genome
interleaveListToGenome l = Genome (Vec.fromList g_g) (Vec.fromList g_ir) l
  where
    (g_g, g_ir) = foldr go ([], []) l
    go l (g, ir)
      | BS.head l == BS.head "g" = ((readG . BS.tail $ l) : g, ir)
      | BS.head l == BS.head "i" = (g, (readG . BS.tail $ l) : ir)
      | otherwise = error patternError
    readG :: (Num a) => BS.ByteString -> a
    readG = fromIntegral . fst . fromJust . BS.readInt

---------------------------------------------------
--      Manipulate Interleave Representation     --
---------------------------------------------------

-- valid paths of a suffix tree must start with a gene
-- must end in a gene
validBeginIR :: BS.ByteString -> Bool
validBeginIR bs = BS.head bs == BS.head "g"

-- valid ByteStrings correspondent to nodes of a suffix tree
-- must end in a gene
validEndIR :: BS.ByteString -> Bool
validEndIR bs = BS.head bs == BS.head "g"

------------------------------------------
--           Operations                 --
------------------------------------------

-- | Replace intergenic regions with zeros
zeroIR :: Genome -> Genome
zeroIR g = Genome vs vi ilr (genomeIsSigned g)
  where
    vs = gstring g
    vi = fmap (const 0) $ irList g
    ilr = Vec.fromList $ interleavelists ls_bs li_bs
    ls_bs = fmap geneToBS . toList $ vs
    li_bs = fmap irToBS . toList $ vi

transposition :: Idx -> Idx -> Idx -> IR -> IR -> IR -> Genome -> Genome
transposition i j k x y z g =
  assert (2 <= i)
    . assert (i < j)
    . assert (j < k)
    . assert (k <= coerce (genomeSize g))
    . assert (0 <= x && x <= (vi ! (coerce i - 2)))
    . assert (0 <= y && y <= (vi ! (coerce j - 2)))
    . assert (0 <= z && z <= (vi ! (coerce k - 2)))
    $ Genome vs' vi' ilr (genomeIsSigned g)
  where
    ilr = Vec.fromList $ interleavelists ls_bs li_bs
    ls_bs = fmap geneToBS . toList $ vs'
    li_bs = fmap irToBS . toList $ vi'
    vs' = Vec.modify updateG vs
    vi' = Vec.modify updateIR vi

    updateG v =
      do
        aux1 <- MVec.clone . MVec.slice (coerce i - 1) (coerce $ j - i) $ v
        aux2 <- MVec.clone . MVec.slice (coerce j - 1) (coerce $ k - j) $ v
        MVec.move (MVec.slice (coerce i - 1) (coerce $ k - j) v) aux2
        MVec.move (MVec.slice (coerce $ i + k - j - 1) (coerce $ j - i) v) aux1

    updateIR v = do
      do
        aux1 <- MVec.clone . MVec.slice (coerce i - 1) (coerce $ j - i) $ v
        aux2 <- MVec.clone . MVec.slice (coerce j - 1) (coerce $ k - j) $ v
        MVec.move (MVec.slice (coerce i - 1) (coerce $ k - j) v) aux2
        MVec.move (MVec.slice (coerce $ i + k - j - 1) (coerce $ j - i) v) aux1
      MVec.write v (coerce i - 1 - 1) (x + y_rest)
      MVec.write v (coerce $ i + k - j - 2) (z + x_rest)
      MVec.write v (coerce k - 1 - 1) (y + z_rest)

    vs = gstring g
    vi = irList g
    x_rest = (vi ! (coerce i - 2)) - x
    y_rest = (vi ! (coerce j - 2)) - y
    z_rest = (vi ! (coerce k - 2)) - z

reversal :: Idx -> Idx -> IR -> IR -> Genome -> Genome
reversal i j x y g =
  assert (2 <= i)
    . assert (i < j)
    . assert (j <= coerce (genomeSize g) - 1)
    . assert (0 <= x && x <= (vi ! (coerce i - 2)))
    . assert (0 <= y && y <= (vi ! (coerce j - 1)))
    $ Genome vs' vi' ilr (genomeIsSigned g)
  where
    ilr = Vec.fromList $ interleavelists ls_bs li_bs
    ls_bs = fmap geneToBS . toList $ vs'
    li_bs = fmap irToBS . toList $ vi'
    vs' = Vec.modify updateG vs
    vi' = Vec.modify updateIR vi

    updateG v = do
      mapM_ (\k -> MVec.swap v (coerce $ i + k - 1) (coerce $ j - k - 1)) [0 .. (j - i + 1) `div` 2 - 1]
      case genomeIsSigned g of
        Unsigned -> pure ()
        Signed -> mapM_ (MVec.modify v invOri . coerce) [i .. j]

    updateIR v = do
      mapM_ (\k -> MVec.swap v (coerce i + coerce k - 1) (coerce j - coerce k - 2)) [0 .. (j - i + 1) `div` 2 - 1]
      MVec.write v (coerce i - 2) (x + y)
      MVec.write v (coerce j - 1) (x_rest + y_rest)

    vs = gstring g
    vi = irList g
    x_rest :: IR
    x_rest = (vi ! (coerce i - 2)) - x
    y_rest :: IR
    y_rest = (vi ! (coerce j - 1)) - y

deletion :: Idx -> Idx -> IR -> Genome -> Genome
deletion i j x g =
  assert (2 <= i)
    . assert (i < j)
    . assert (j <= coerce (genomeSize g))
    . assert (0 <= x && x <= (vi ! (coerce i - 2)) + (vi ! (coerce j - 2)))
    $ Genome vs' vi' ilr (genomeIsSigned g)
  where
    ilr = Vec.fromList $ interleavelists ls_bs li_bs
    ls_bs = fmap geneToBS . toList $ vs'
    li_bs = fmap irToBS . toList $ vi'
    vs' =
      Vec.slice 0 (coerce i - 1) vs
        Vec.++ Vec.slice (coerce j - 1) (Vec.length vs - coerce j + 1) vs
    vi' =
      Vec.slice 0 (coerce i - 2) vi
        Vec.++ Vec.fromList [x]
        Vec.++ if coerce j == genomeSize g
          then Vec.empty
          else Vec.slice (coerce j - 1) (Vec.length vi - coerce j + 1) vi

    vs = gstring g
    vi = irList g

-- Genome to be inserted must be open (n genes and n+1 itergenic regions)
insertion :: Idx -> Genome -> Genome -> Genome
insertion i x g =
  assert (1 <= i)
    . assert (i <= coerce (genomeSize g) - 1)
    $ Genome vs' vi' ilr (genomeIsSigned g)
  where
    ilr = Vec.fromList $ interleavelists ls_bs li_bs
    ls_bs = fmap geneToBS . toList $ vs'
    li_bs = fmap irToBS . toList $ vi'
    vs' =
      Vec.slice 0 (coerce i) vs
        Vec.++ vsx
        Vec.++ Vec.slice (coerce i) (Vec.length vs - coerce i) vs
    vi' =
      ( if i == 1
          then Vec.empty
          else Vec.slice 0 (coerce i - 1) vi
      )
        Vec.++ vix
        Vec.++ Vec.slice (coerce i) (Vec.length vi - coerce i) vi

    vs = gstring g
    vi = irList g
    vsx = gstring x
    vix = irList x

combineGenomes :: IR -> Genome -> Genome -> Genome
combineGenomes ir g h = Genome (gstring g <> gstring h) (irList g <> Vec.cons ir (irList h)) (Vec.concat [interleaveListRepresentation g, Vec.singleton (irToBS ir), interleaveListRepresentation h]) (genomeIsSigned g <> genomeIsSigned h)

combineGenomesL :: [IR] -> [Genome] -> Genome
combineGenomesL irs gs = Genome vs vi ilr (genomeIsSigned $ head gs)
      where
          vs = Vec.concat (map gstring gs)
          vi = Vec.concat $ interleavelists (map irList gs) (map pure irs)
          ilr = Vec.fromList $ interleavelists ls_bs li_bs
          ls_bs = fmap geneToBS . toList $ vs
          li_bs = fmap irToBS . toList $ vi

breakGenome :: Genome -> Idx -> (IR, Genome, Genome)
breakGenome g idx = (irList g ! (coerce idx - 1), sliceGenome g 1 idx, sliceGenome g (idx + 1) (coerce $ genomeSize g))

sliceGenome :: Genome -> Idx -> Idx -> Genome
sliceGenome g i j = Genome vs' vi' ilr' (genomeIsSigned g)
  where
    ilr' = Vec.slice (coerce $ 2 * (i -1)) (coerce $ 2 * (j - i) + 1) $ interleaveListRepresentation g
    vs' = Vec.slice (coerce $ i - 1) (coerce $ j - i + 1) $ gstring g
    vi' = Vec.slice (coerce $ i - 1) (coerce $ j - i) $ irList g

------------------------------------------
--           Descriptions               --
------------------------------------------

irByIdx :: Genome -> Idx -> IR
irByIdx g idx = irList g ! (coerce idx - 1)

duoByIdx :: Genome -> Idx -> Duo
duoByIdx g idx = Duo (a1, a2) ir idx (genomeSize g) (genomeIsSigned g)
  where
    i = coerce idx
    vs = gstring g
    vi = irList g
    a1 = vs ! (i -1)
    a2 = vs ! i
    ir = vi ! (i -1)

occurence :: Genome -> Gene -> Int
occurence g a = sum . fmap (\x -> if x == a then 1 else 0) . gstring $ g

occurenceMap :: Genome -> IntMap Int
occurenceMap = IntMap.fromListWith (+) . (`zip` [1,1..]) . map abs . coerce . toList . gstring

occurenceMapLookup :: Gene -> IntMap Int -> Maybe Int
occurenceMapLookup a = IntMap.lookup (abs . coerce $ a)

occurenceMapAdjust :: (Int -> Int) -> Gene -> IntMap Int -> IntMap Int
occurenceMapAdjust f a = IntMap.adjust f (abs . coerce $ a)

occurenceMax :: Genome -> Int
occurenceMax = maximum . occurenceMap

alphabet :: Genome -> [Gene]
alphabet = unique . toList . gstring

balanced :: Genome -> Genome -> Bool
balanced g h = balancedGenes g h && balancedIR g h
    where
        balancedGenes g h =
          (List.sort . map canonicOri . Vec.toList $ gstring g)
            == (List.sort . map canonicOri . Vec.toList $ gstring h)
        balancedIR g h = sum (irList g) == sum (irList h)

genomeSize :: Genome -> Size
genomeSize = Size . Vec.length . gstring

-- Maximum value of a gene
geneMaxValue :: Genome -> Gene
geneMaxValue = maximum . fmap (\v -> if v == maxBound then 0 else abs v) . gstring

------------------------------------------
--           Sub-Genome                 --
------------------------------------------

subGenome :: Genome -> Genome -> Bool
subGenome x g = genomeIsSigned x == genomeIsSigned g && or (subGenomeFind x g)

-- | For each gene of a genome g, return if an occurence of
--  a nonempty genome x starts in that gene
subGenomeFind :: Genome -> Genome -> Vector Bool
subGenomeFind x g =
  fmap (\i -> (xs `isPrefixOfV` Vec.drop i gs) && (xi `isPrefixOfV` Vec.drop i gi)) inds
  where
    inds = Vec.elemIndices (Vec.head $ gstring x) gs
    gs = gstring g
    gi = irList g
    xs = gstring x
    xi = irList x
    isPrefixOfV v1 v2 = length v1 <= length v2 && v1 == Vec.unsafeTake (length v1) v2

allSubGenomes :: Genome -> [Genome]
allSubGenomes g = go 0 1 []
  where
    vs = gstring g
    vi = irList g
    n = coerce $ genomeSize g
    go i k acc
      | i > n = acc
      | i + k > n = go (i + 1) 1 acc
      | otherwise =
        let g' = Genome (Vec.slice i k vs) (Vec.slice i (k -1) vi) (Vec.slice (2 * i) (2 * k - 1) $ interleaveListRepresentation g) (genomeIsSigned g)
         in go i (k + 1) (g' : acc)

maybeEstimateITD :: Genome -> Genome -> Maybe Dist
maybeEstimateITD g h =
  if balanced g h
    then Just $ estimateITD g h
    else Nothing

-- | Count the number of occurrences of a nonempty genome x
--  in a genome g
subGenCount :: Genome -> Genome -> Int
subGenCount x g = sum . fmap (\b -> if b then 1 else 0) $ subGenomeFind x g

-- weigth :: Genome -> Genome -> Genome -> Int
-- weigth g h x = subGenCount x g - subGenCount x h

------------------------------------------
--              Distances               --
------------------------------------------

estimateITD :: Genome -> Genome -> Dist
estimateITD g h = undefined

----------------------------------
--        Genome Map            --
----------------------------------

data MapType = Direct | Reverse

data GenomeMap v = GM MapType (IntMap [(Genome, v)])

instance Show v => Show (GenomeMap v) where
  show (GM _ m) = show . concat . IntMap.elems $ m

newtype GeneListForMap = MkGeneListForMap {unGeneListForMap :: [Gene]} deriving newtype (Eq, Show)

instance Orientable GeneListForMap where
  getOri = getOri . head . unGeneListForMap
  invOri = MkGeneListForMap . map invOri . unGeneListForMap

gmEmpty :: GenomeMap v
gmEmpty = GM Direct IntMap.empty

gmEmptyRev :: GenomeMap v
gmEmptyRev = GM Reverse IntMap.empty

gmLookup :: (Orientable v) => Genome -> GenomeMap v -> Maybe v
gmLookup g m@(GM Direct _) = gmLookup_ g k m
  where
    k = coerce . Vec.head . gstring $ g
gmLookup g m@(GM Reverse _) =
  case getOri g of
    LR -> gmLookup_ g' k m
    RL -> invOri <$> gmLookup_ g' k m
  where
    g' = canonicOri g
    k = coerce . Vec.head . gstring $ g'

gmAlter :: forall v. Genome -> (Maybe v -> v) -> GenomeMap v -> GenomeMap v
gmAlter g f (GM mtype m) = GM mtype (IntMap.alter f' k m)
  where
    f' :: Maybe [(Genome,v)] -> Maybe [(Genome,v)]
    f' Nothing = Just [(g', f Nothing)]
    f' (Just ls) = f'' [] ls
    f'' :: [(Genome, v)] -> [(Genome, v)] -> Maybe [(Genome, v)]
    f'' acc [] = Just $ reverse acc ++ [(g', f Nothing)]
    f'' acc ((g',val):ls) =
        if g' == g
           then Just $ reverse acc ++ (g', f (Just val)) : ls
           else f'' ((g',val):acc) ls
    g' = case mtype of Direct -> g; Reverse -> canonicOri g
    k = coerce . Vec.head . gstring $ g'

gmLookup_ :: Genome -> Int -> GenomeMap v -> Maybe v
gmLookup_ g k (GM _ m) = fmap snd $ IntMap.lookup k m >>= List.find ((==) g . fst)

gmInsert :: Genome -> v -> GenomeMap v -> GenomeMap v
gmInsert g v m@(GM mtype _) = gmInsert_ g' k v m
  where
    g' = case mtype of Direct -> g; Reverse -> canonicOri g
    k = coerce . Vec.head . gstring $ g'

gmInsert_ :: Genome -> Int -> v -> GenomeMap v -> GenomeMap v
gmInsert_ g k val (GM t m) = GM t m'
  where
    m' = case IntMap.lookup k m of
      Nothing -> IntMap.insert k [(g, val)] m
      Just _ -> IntMap.adjust ([(g, val)] ++) k m

-- if not found (return False and insert) else (return True)
gmLookupInsert :: (Orientable v) => Genome -> v -> GenomeMap v -> (GenomeMap v, v, Bool)
gmLookupInsert g val m@(GM Direct _) =
  case gmLookup_ g k m of
    Nothing -> (gmInsert_ g k val m, val, False)
    Just v -> (m, v, True)
  where
    k = coerce . Vec.head . gstring $ g
gmLookupInsert g val m@(GM Reverse _) =
  case gmLookup_ g' k m of
    Nothing -> (gmInsert_ g' k val m, val, False)
    Just v -> (m, case getOri g of LR -> v; RL -> invOri v, True)
  where
    g' = canonicOri g
    k = coerce . Vec.head . gstring . canonicOri $ g'
