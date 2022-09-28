{-# LANGUAGE TemplateHaskell #-}

module BalancerCheck (tests) where

import Balancer
import LocalBase
import Genomes (occurenceMap)
import GenomesCheck (genGenomeWithSign, genSign)
import qualified Data.IntMap as IntMap
import qualified Data.List as List
import Hedgehog

prop_resultIsBalanced :: Property
prop_resultIsBalanced = property $ do
  signed <- forAll genSign
  g <- forAll $ genGenomeWithSign signed
  h <- forAll $ genGenomeWithSign signed
  let (g', h') = reduceGenesForDeletion g h
      occMap = IntMap.unionWith (+) (occurenceMap g') (negate <$> occurenceMap h')
      (bal,unbal) = List.partition (\(k,v) -> v == 0 && k /= maxBound) $ IntMap.assocs occMap
  assert (null bal || null unbal || maximum bal < minimum unbal)

tests :: IO Bool
tests = checkSequential $$(discover)
