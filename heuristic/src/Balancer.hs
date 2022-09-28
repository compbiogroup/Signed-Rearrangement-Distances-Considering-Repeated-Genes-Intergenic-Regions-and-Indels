module Balancer where

import Control.Exception (assert)
import Data.Coerce (coerce)
import Data.Foldable (toList)
import Data.IntMap.Strict (IntMap)
import qualified Data.IntMap.Strict as IntMap
import Data.IntSet (IntSet)
import qualified Data.IntSet as IntSet
import Data.List (foldl', mapAccumL)
import Data.Sequence (Seq ((:<|), (:|>)))
import qualified Data.Sequence as Seq
import Genomes as G
import LocalBase

data Interval = Interval
  { b_start :: Idx,
    b_end :: Idx,
    b_size :: Size,
    b_occMap :: IntMap Int
  } deriving (Show)

-- | Replace some blocks with new genes. Those genes will mark the regions that must be deleted.
-- the resulting genomes are balanced without the marks.
reduceGenesForDeletion :: Genome -> Genome -> (Genome, Genome)
reduceGenesForDeletion g h = (g', h')
  where
    (deleterG, deleterH) = markGenesForDelition g h
    gmax = max (geneMaxValue g) (geneMaxValue h) + 1
    (sign, g_lg, g_li) = toLists False g
    (_, h_lg, h_li) = toLists False h
    (g_lg', g_li', gmax') = copyGene deleterG ([], [], gmax, 1) g_lg g_li
    g' = fromLists False (genomeIsSigned g) g_lg' g_li'
    (h_lg', h_li', _) = copyGene deleterH ([], [], gmax', 1) h_lg h_li
    h' = fromLists False (genomeIsSigned h) h_lg' h_li'

    copyGene :: IntSet -> ([Gene], [IR], Gene, Idx) -> [Gene] -> [IR] -> ([Gene], [IR], Gene)
    copyGene deleter = copyGene_
      where
        copyGene_ _ [] [] = error "Error in reduceGenesForDeletion, gstring and ir list cannot have same size."
        copyGene_ _ [] (_ : _) = error "Error in reduceGenesForDeletion, ir_list can not be bigger than gstring."
        copyGene_ (acc_lg, acc_li, val, i) genes@(a : as) [] = assert (not (coerce i `IntSet.member` deleter)) (reverse (a : acc_lg), reverse acc_li, val)
        copyGene_ (acc_lg, acc_li, val, i) genes@(a : as) ir_list@(ir : irs) = if coerce i `IntSet.member` deleter then deleteGene (acc_lg, acc_li, val, ir, i) genes ir_list else copyGene_ (a : acc_lg, ir : acc_li, val, i + 1) as irs

        deleteGene :: ([Gene], [IR], Gene, IR, Idx) -> [Gene] -> [IR] -> ([Gene], [IR], Gene)
        deleteGene _ [] [] = error "Error in reduceGenesForDeletion, gstring and ir list cannot have same size (del)."
        deleteGene _ [] (_ : _) = error "Error in reduceGenesForDeletion, ir_list can not be bigger than gstring (del)."
        deleteGene (acc_lg, acc_li, val, old_ir, i) genes@(a : as) [] = copyGene_ (val : acc_lg, old_ir : acc_li, val + 1, i) genes []
        deleteGene (acc_lg, acc_li, val, old_ir, i) genes@(a : as) ir_list@(ir : irs) = if coerce i `IntSet.member` deleter then deleteGene (acc_lg, acc_li, val, ir, i + 1) as irs else copyGene_ (val : acc_lg, old_ir : acc_li, val + 1, i) genes ir_list

    markGenesForDelition :: Genome -> Genome -> (IntSet, IntSet)
    markGenesForDelition g h = (deleterG, deleterH)
      where
        (deleterG, occMapG) = markIntervals IntSet.empty occMap0 1 g_lg
        (deleterH, _) = markIntervals IntSet.empty occMapG (-1) h_lg
        occMap0 = IntMap.unionWith (+) (occurenceMap g) (negate <$> occurenceMap h)
        -- Check each interval and select the bigger one,
        -- the tiebraker is the position of the interval (leftmost selected).
        markIntervals :: IntSet -> IntMap Int -> Int -> [Gene] -> (IntSet, IntMap Int)
        markIntervals deleter occMap negOrPos lg =
          if i_start == 0
            then (deleter, occMap)
            else markIntervals deleter' occMap' negOrPos lg
          where
            deleter' = foldl' (flip IntSet.insert) deleter (coerce [i_start .. i_end] :: [Int])
            (Interval i_start i_end i_size occMap', _, _) = foldl' checkGeneForDelition (Interval 0 0 0 occMap, Interval 0 0 0 occMap, False) $ zip [1..] lg

            checkGeneForDelition :: (Interval, Interval, Bool) -> (Idx,Gene) -> (Interval, Interval, Bool)
            checkGeneForDelition (interval_max, interval_current, inB) (pos,a) = (interval_max', interval_current', inB')
              where
                (interval_current', inB') =
                  case occurenceMapLookup a (b_occMap interval_current) of
                    Nothing -> error "Error in reduceGenesForDeletion, key must be in map"
                    Just occ ->
                      let occMap_updated = occurenceMapAdjust (subtract negOrPos) a (b_occMap interval_current)
                       in if coerce pos `IntSet.notMember` deleter && negOrPos * occ > 0
                            then if inB
                                    then (Interval (b_start interval_current) (b_end interval_current + 1) (b_size interval_current + 1) occMap_updated, True)
                                    else (Interval pos pos 1 occMap_updated, True)
                            else (Interval 0 0 0 occMap, False)
                interval_max' =
                  if inB' && b_size interval_max < b_size interval_current'
                    then interval_current'
                    else interval_max
