{-# LANGUAGE TemplateHaskell #-}

module GenomesCheck (tests, genIR, genGenome, genGenomeWithSign, genSign, getSubGenome, GenomesCheck.rearrangeGenome) where

import Control.Monad.Random (evalRandIO)
import Control.Monad.Trans.Class (lift)
import Control.Monad.Trans.Maybe (MaybeT (MaybeT), runMaybeT)
import Data.Coerce (coerce)
import qualified Data.List as List
import Data.Maybe (fromJust, isJust)
import qualified Data.Vector as Vec
import Genomes
import Hedgehog
import qualified Hedgehog.Gen as Gen
import qualified Hedgehog.Range as Range
import LocalBase

genGene :: Gen Gene
genGene = intToGene <$> Gen.int (Range.linear 1 100)

genIR :: Gen IR
genIR = coerce <$> Gen.int (Range.linear 0 100)

genSign :: Gen Sign
genSign = Gen.enum Signed Unsigned

-- | Generate an empty genome
genGenome :: Gen Genome
genGenome = do
  signed <- genSign
  genGenomeWithSign signed

genGenomeWithSign :: Sign -> Gen Genome
genGenomeWithSign signed = do
  size <- Gen.int (Range.linear 1 100)
  coins <- Gen.list (Range.singleton $ size - 1) Gen.bool
  ls <-
    ( case signed of
        Unsigned -> id
        Signed -> zipWith swaps coins
      )
      <$> Gen.list (Range.singleton $ size - 1) genGene
  li <- Gen.list (Range.singleton size) genIR
  return $ fromLists True signed ls li
  where
    swaps b v = if b then v else invOri v

rearrangeGenome :: Genome -> Gen Genome
rearrangeGenome g = do
  let (_, ls, li) = toLists True g
      s = sum li
  ls' <- Gen.shuffle ls
  x <-
    fmap List.sort . Gen.list (Range.singleton (coerce $ genomeSize g - 2))
      . fmap IR
      . Gen.int
      $ Range.constant 0 (coerce s)
  let li' = zipWith (-) (x ++ [s]) (0 : x)
  return $ fromLists True (genomeIsSigned g) ls' li'

getSubGenome :: Genome -> Gen Genome
getSubGenome g = fmap fromJust . Gen.filter isJust . runMaybeT $ do
  idx <- MaybeT $ getIdx g
  let (_, gl, gr) = breakGenome g idx
  g' <- lift $ Gen.element [gl, gr]
  idx' <- MaybeT $ getIdx g'
  let (_, gl', gr') = breakGenome g' idx'
  lift $ Gen.element [gl', gr']

getIdx :: Genome -> Gen (Maybe Idx)
getIdx g =
  if genomeSize g == 1
    then pure Nothing
    else fmap (Just . coerce) . Gen.int $ Range.constant 1 (coerce $ genomeSize g - 1)

prop_balancedAfterRearrange :: Property
prop_balancedAfterRearrange =
  property $ do
    g <- forAll genGenome
    h <- evalIO . evalRandIO $ Genomes.rearrangeGenome g
    assert $ balanced g h

prop_combineIsAssociative :: Property
prop_combineIsAssociative =
  property $ do
    g <- forAll genGenome
    h <- forAll genGenome
    k <- forAll genGenome
    ri1 <- forAll genIR
    ri2 <- forAll genIR
    combineGenomes ri1 g (combineGenomes ri2 h k) === combineGenomes ri2 (combineGenomes ri1 g h) k

prop_breakVsCombine :: Property
prop_breakVsCombine =
  property $ do
    g <- forAll . Gen.filter (\x -> genomeSize x >= 2) $ genGenome
    break <- forAll . fmap fromJust $ getIdx g
    (ir, gl, gr) <- eval $ breakGenome g break
    g === combineGenomes ir gl gr

prop_subgenomeAndEquality :: Property
prop_subgenomeAndEquality =
  property $ do
    g <- forAll genGenome
    h <- forAll genGenome
    (subGenome g h && subGenome h g) === (g == h)

prop_subgenomeIsTransitive :: Property
prop_subgenomeIsTransitive =
  property $ do
    g <- forAll genGenome
    h <- forAll genGenome
    k <- forAll genGenome
    if subGenome g h && subGenome h k
      then assert (subGenome g k)
      else success

prop_subGenCountIsCorrect :: Property
prop_subGenCountIsCorrect =
  property $ do
    g <- forAll genGenome
    x <- forAll $ getSubGenome g
    l <- forAll . pure $ allSubGenomes g
    subGenCount x g === length [y | y <- l, y == x]

prop_toFromBSIso :: Property
prop_toFromBSIso =
  property $ do
    g <- forAll genGenome
    g === (uncurry (readGenome False (genomeIsSigned g)) . writeGenome False $ g)
    g === (uncurry (readGenome True (genomeIsSigned g)) . writeGenome True $ g)

prop_toFromListIso :: Property
prop_toFromListIso =
  property $ do
    g <- forAll genGenome
    g === ((\(_, x, y) -> fromLists False (genomeIsSigned g) x y) . toLists False $ g)
    g === ((\(_, x, y) -> fromLists True (genomeIsSigned g) x y) . toLists True $ g)

prop_interleaveListIso :: Property
prop_interleaveListIso =
  property $ do
    g <- forAll genGenome
    g === (flip interleaveListToGenome (genomeIsSigned g) . interleaveListRepresentation $ g)

tests :: IO Bool
tests = checkSequential $$(discover)
