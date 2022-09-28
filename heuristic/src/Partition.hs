{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}

-- |
-- Module      : Partition
-- Description : Representation and construction of genome partitions
-- Copyright   : (c) Gabriel Siqueira, 2021
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
--
-- A partition of two genomes G and H is essentially composed of two genome sequences S and P, such that:
--     - the sequence S when combined gives us G
--     - the sequence P when combined gives us H
--     - it is possible to rearrange S to obtain P
module Partition
  ( Partition,
    PartitionType (..),
    getPartition,
    validPartition,
    breakpoints,
    blocks,
    mapToPerm,
    reduced,
    cost,
    -- Partition.weigth,
    sizeTmin,
  )
where

import Control.Arrow ((***))
import Control.Exception (assert)
import Data.ByteString (ByteString)
import qualified Data.ByteString.Char8 as BS
import Data.Coerce (coerce)
import Data.Foldable (toList)
import Data.HashSet (HashSet)
import qualified Data.HashSet as HashSet
import qualified Data.List as List
import qualified Data.Map as Map
import Data.Maybe
  ( catMaybes,
    fromJust,
    fromMaybe,
    listToMaybe,
    mapMaybe,
  )
import Data.Sequence (Seq ((:<|), (:|>)))
import qualified Data.Sequence as Seq
import qualified Data.Set as Set
import Data.Vector (Vector, (!))
import qualified Data.Vector as Vec
import Debug.Trace
import GTree
import Genomes as G
import LocalBase

-- | Representation of a partition
newtype Partition = ValidPartition PartialPartition deriving newtype (Show)

data PartialPartition = PPartition {partType :: PartitionType, gseq :: Seq Genome, hseq :: Seq Genome, gbps :: Seq IR, hbps :: Seq IR} deriving (Show)

-- | A partial partition may be invalid
makePartialPartition :: PartitionType -> Genome -> Genome -> PartialPartition
makePartialPartition ptype g h = PPartition ptype (pure g) (pure h) mempty mempty

-- | Recives: a partial partition and an indication of the original genome and
-- index of an occurrence of a duo containing the breakpoint
--   Returns: a partial partition with the new breakpoint
addBreakpoint :: PartialPartition -> GenomePosition -> PartialPartition
addBreakpoint part duoPos =
  case duoPos of
    (G idx _ _) ->
      let (gseq', gbps') = aux (gseq part) (gbps part) idx
       in part {gseq = gseq', gbps = gbps'}
    (H idx _ _) ->
      let (hseq', hbps') = aux (hseq part) (hbps part) idx
       in part {hseq = hseq', hbps = hbps'}
  where
    aux :: Seq Genome -> Seq IR -> Idx -> (Seq Genome, Seq IR)
    aux kseq kbps idx = (kseq', kbps')
      where
        (kseq_front, x, kseq_back, kbps_front, kbps_back, idx') = go Seq.Empty kseq Seq.Empty kbps idx

        (ir, xl, xr) = breakGenome x idx'
        kseq' = if idx' == 0 then kseq else kseq_front <> (xl :<| xr :<| kseq_back)
        kbps' = if idx' == 0 then kbps else kbps_front <> (ir :<| kbps_back)

        -- Go through the elements of the genome sequence until it finds the genome that must
        -- be broken. Returns the first part of the sequence, the genome to be broken, the last
        -- part of the sequence, the two parts of the intergenic region sequence, and an index
        -- indicating where the genome must be broken.
        go ::
          Seq Genome ->
          Seq Genome ->
          Seq IR ->
          Seq IR ->
          Idx ->
          (Seq Genome, Genome, Seq Genome, Seq IR, Seq IR, Idx)
        go ys (x :<| xs) ybs xbss@(xb :<| xbs) idx' =
          let size = coerce (genomeSize x)
           in if idx' >= size
                then go (ys :|> x) xs (ybs :|> xb) xbs (idx' - size)
                else (ys, x, xs, ybs, xbss, idx')
        go ys (x :<| xs) ybs xbss idx' = (ys, x, xs, ybs, xbss, idx')
        go _ Seq.Empty _ _ _ = error patternError

-- | 2k-approximation for the intergenic partition problem
getPartition :: PartitionType -> Genome -> Genome -> Partition
getPartition ptype g h = ValidPartition $ getPartition_ (makePartialPartition ptype g h) update ptype g h
  where
    update part _ breaksPos = foldl addBreakpoint part (map snd breaksPos)

listGenomes :: PartitionType -> Genome -> Genome -> [Genome]
listGenomes = getPartition_ [] update
  where
    update l x _ = x : l

-- | Produce a value of type `a` by selecting elements from a Tmin set.
-- The update function must receive a value of type 'a', the selected genome
-- and a list of breakpoints to be inserted in the partition.
getPartition_ :: forall a. a -> (a -> Genome -> [(Maybe Genome, GenomePosition)] -> a) -> PartitionType -> Genome -> Genome -> a
getPartition_ acc0 updateAcc ptype g h = go acc0 tmin0 breaks0
  where
    tmin0 = makeTmin ptype g h
    breaks0 = Breaks ptype HashSet.empty
    go :: a -> Tmin -> Breaks -> a
    go acc tmin breaks =
      case getGenome tmin of
        Nothing -> acc
        Just (x, gp) -> assert (tmin' /= tmin) $ go acc' tmin' breaks'
          where
            tmin' =
              case ptype of
                MCISP -> tmin1
                RMCISP -> tmin2
            (breaksPos, breaks') = getBreak g h gp breaks x
            tmin1 = foldr (updateTmin g h) tmin breaksPos
            tmin2 =
              foldr
                ( updateTmin (invOri g) (invOri h)
                    . (invOri *** invOri)
                )
                tmin1
                breaksPos
            acc' = updateAcc acc x breaksPos

-- | A partition (s,p) is valid if
--  it is possible to rearrange s to obtain p.
validPartition :: Partition -> Bool
validPartition = uncurry balanced . reduced

breakpoints :: Partition -> (Seq IR, Seq IR)
breakpoints (ValidPartition part) = (gbps part, hbps part)

blocks :: Partition -> (Seq Genome, Seq Genome)
blocks (ValidPartition part) = (gseq part, hseq part)

-- | Converts a partition into a pair of genomes representing a mapping into
-- a permutation compatible with the reduced genomes correspondent to the partition
mapToPerm :: Genome -> Genome -> Partition -> (Genome, Genome)
mapToPerm g_bal h_bal (ValidPartition part) = (gr, hr)
  where
    (_,_,g_bal_ir) = toLists False g_bal
    (_,_,h_bal_ir) = toLists False h_bal
    gr = fromLists False sign gr_ls g_bal_ir
    hr = fromLists False sign hr_ls h_bal_ir
    (gr_ls, hr_ls) = (genomesToUniqueList ggs, genomesToUniqueList hhs)
    ggs = gseq part
    hhs = hseq part
    (sign, gmEmptyX) = case partType part of
      MCISP -> (Unsigned, gmEmpty)
      RMCISP -> (Signed, gmEmptyRev)

    -- Produce list of unique characters correspondent to genes
    genomesToUniqueList :: Seq Genome -> [Gene]
    genomesToUniqueList = genomesToUniqueList' [] m_orig . toList
    genomesToUniqueList' :: [Gene] -> GenomeMap GeneListForMap -> [Genome] -> [Gene]
    genomesToUniqueList' acc m [] = reverse acc
    genomesToUniqueList' acc m (g:gs) = genomesToUniqueList' ([v..v + n - 1] ++ acc) m' gs
        where 
            (_,lg,_) = toLists False g
            n = intToGene . coerce . genomeSize $ g
            v = head . unGeneListForMap . fromMaybe (error "Error on mapToPerm (1).") $ gmLookup g m
            m' = gmAlter g (\old -> case old of
                                      Nothing -> error "Error on mapToPerm (2)."
                                      Just (MkGeneListForMap (x:xs)) -> MkGeneListForMap xs) m
    m_orig = fst $ foldl addGenome (gmEmptyX, 1 :: Gene) (hhs Seq.>< ggs)
    addGenome :: (GenomeMap GeneListForMap, Gene) -> Genome -> (GenomeMap GeneListForMap, Gene)
    addGenome (m, count) g = (m', count')
      where
        n = intToGene . coerce . genomeSize $ g
        count' = count + n
        m' = gmAlter g (\old -> case old of
                                  Nothing -> MkGeneListForMap [count]
                                  Just (MkGeneListForMap l) -> MkGeneListForMap (count:l)) m

-- | Converts a partition into a pair of genomes representing the reduced genomes
-- correspondent to the partition
reduced :: Partition -> (Genome, Genome)
reduced (ValidPartition part) = (gr, hr)
  where
    gr = fromLists False sign gr_ls (toList $ gbps part)
    hr = fromLists False sign hr_ls (toList $ hbps part)
    (gr_ls, hr_ls) = (genomesToGenes ggs, genomesToGenes hhs)
    ggs = gseq part
    hhs = hseq part
    (sign, gmEmptyX) = case partType part of
      MCISP -> (Unsigned, gmEmpty)
      RMCISP -> (Signed, gmEmptyRev)

    genomesToGenes = map genomeToGene . toList
    genomeToGene g = fromMaybe (error "Error on reduced.") $ gmLookup g m
    m = fst $ foldl addGenome (gmEmptyX, 1 :: Gene) (hhs Seq.>< ggs)
    addGenome (m, count) g = (m', count')
      where
        count' = if keepOld then count else count + 1
        (m', _, keepOld) = gmLookupInsert g count m

-- | The cost of a partition (s,p) is the number of
--  breakpoints from s
cost :: Partition -> Int
cost = length . fst . breakpoints

-- weigth :: Partition -> Genome -> Int
-- weigth part x = sumCounts bls1 - sumCounts bls2
--   where
--     sumCounts = sum . fmap (subGenCount' x)
--     subGenCount' x g =
--         case partType part of
--           MCISP -> subGenCount x g
--           RMCISP -> subGenCount x g + subGenCount (invOri x) g
--     (bls1, bls2) = blocks part

sizeTmin :: PartitionType -> Genome -> Genome -> Size
sizeTmin ptype g h = Size . length $ listGenomes ptype g h

-- | Breaks represents a set of breakpoints from a partial partition
data Breaks = Breaks PartitionType (HashSet (Gene, Gene)) deriving (Show)

isBreak :: Duo -> Breaks -> Bool
isBreak duo (Breaks partType brks) = (genePair . canonicOri $ duo) `HashSet.member` brks

-- Position of the break that must be inserted, if the genome has size one two breaks are returned (breaks around the genome).
getBreak :: Genome -> Genome -> GenomePosition -> Breaks -> Genome -> ([(Maybe Genome, GenomePosition)], Breaks)
getBreak g h gp brks@(Breaks partType set) x =
  if genomeSize x == 1
    then case gp of
      (G gidx _ ori) ->
        let breaks_list = (\idx -> G idx (genomeSize g - 2) ori) <$> [gidx - 1, gidx]
            duo1 = genePair . canonicOri $ duoByIdx g (gidx - 1)
            duo2 = genePair . canonicOri $ duoByIdx g gidx
         in ( zip [Nothing, Just x] breaks_list,
              Breaks partType (HashSet.insert duo2 (HashSet.insert duo1 set))
            )
      (H gidx _ ori) ->
        let breaks_list = (\idx -> H idx (genomeSize h - 2) ori) <$> [gidx - 1, gidx]
            duo1 = genePair . canonicOri $ duoByIdx h (gidx - 1)
            duo2 = genePair . canonicOri $ duoByIdx h gidx
         in ( zip [Nothing, Just x] breaks_list,
              Breaks partType (HashSet.insert duo2 (HashSet.insert duo1 set))
            )
    else (,brks') $ case gp of
      (G gidx _ ori) -> [(Nothing, G (gidx + didx - 1) (genomeSize g - 2) ori)]
      (H gidx _ ori) -> [(Nothing, H (gidx + didx - 1) (genomeSize h - 2) ori)]
  where
    didx = duoIdx break
    (break, brks') = case dropWhile (not . (`isBreak` brks)) ls of
      [] -> (head ls, Breaks partType (HashSet.insert (genePair . canonicOri $ head ls) set))
      (l : _) -> (l, brks)
    ls = duosList x
