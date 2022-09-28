{-# LANGUAGE TemplateHaskell #-}

module PartitionCheck (tests) where

import Control.Concurrent (threadDelay)
import Control.Concurrent.Async.Lifted (race)
import Control.Monad.IO.Class (liftIO)
import Data.Coerce (coerce)
import Data.Foldable (toList)
import Genomes as G hiding (rearrangeGenome)
import GenomesCheck (genGenome, genGenomeWithSign, genSign, rearrangeGenome)
import Hedgehog
import qualified Hedgehog.Gen as Gen
import qualified Hedgehog.Range as Range
import LocalBase
import Partition as P

withTimeLimit :: Int -> TestT IO a -> TestT IO a
withTimeLimit timeout v = do
  result <-
    race
      (liftIO $ threadDelay timeout)
      v
  case result of
    Left () -> fail "Timeout exceeded"
    Right x -> pure x

genPartitionType :: Gen PartitionType
genPartitionType = Gen.enumBounded

genGenomes :: Gen (Genome, Genome, PartitionType)
genGenomes = do
  signed <- genSign
  g <- genGenomeWithSign signed
  h <- genGenomeWithSign signed
  pt <- genPartitionType
  return (g,h,pt)

genBalancedGenomes :: Gen (Genome, Genome, PartitionType)
genBalancedGenomes = do
  g <- genGenome
  h <- rearrangeGenome g
  pt <- genPartitionType
  return (g,h,pt)

genBalancedPartition :: Gen Partition
genBalancedPartition = do
  (g,h,pt) <- genBalancedGenomes
  return $ getPartition pt g h

prop_partitionSequences :: Property
prop_partitionSequences =
  property $ do
    (g,h,pt) <- forAll genGenomes
    part <- test . withTimeLimit 1000000 . eval $ getPartition pt g h
    let (bps1, bps2) = breakpoints part
        (bls1, bls2) = blocks part
    combineGenomesL (toList bps1) (toList bls1) === g
    combineGenomesL (toList bps2) (toList bls2) === h

prop_partitionIsValid :: Property
prop_partitionIsValid =
  property $ do
    part <- forAll genBalancedPartition
    assert (validPartition part)

prop_costIsSimetric :: Property
prop_costIsSimetric =
  property $ do
    (g,h,pt) <- forAll genGenomes
    let partGH = getPartition pt g h
    let partHG = getPartition pt g h
    cost partGH === cost partHG

prop_costIsCorrect :: Property
prop_costIsCorrect =
  property $ do
    part <- forAll genBalancedPartition
    cost part === (subtract 1 . length . fst . blocks $ part)

-- prop_makeTminIsCorrect :: Property
-- prop_makeTminIsCorrect =
--     property $ do
--         return undefined

prop_lowerBound :: Property
prop_lowerBound =
  property $ do
    (g,h,pt) <- forAll genBalancedGenomes
    let part = getPartition pt g h
    assert $ coerce (sizeTmin pt g h) <= 2 * cost part

prop_upperBound :: Property
prop_upperBound =
  property $ do
    (g,h,pt) <- forAll genBalancedGenomes
    let part = getPartition pt g h
        k = occurenceMax g
    assert $ 2 * k * coerce (sizeTmin pt g h) >= 2 * cost part

tests :: IO Bool
tests = checkSequential $$(discover)
