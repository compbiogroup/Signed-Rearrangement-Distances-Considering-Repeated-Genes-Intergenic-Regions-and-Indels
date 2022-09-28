{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TupleSections #-}

-- |
-- Module      : GTree
-- Description : Representation and construction of a sub-genome tree
-- Copyright   : (c) Gabriel Siqueira, 2021
-- License     : BSD3
-- Maintainer  : gabriel.gabrielhs@gmail.com
--
-- This module implements a suffix tree based structure to represent common sub-genomes of a pair of genomes.
module GTree
  ( Tmin,
    GenomePosition (..),
    PartitionType (..),
    makeTmin,
    getGenome,
    updateTmin,
  )
where

import Control.Arrow (second)
import Control.Exception (assert)
import Data.ByteString (ByteString)
import qualified Data.ByteString.Char8 as BS
import Data.Coerce (coerce)
import Data.Foldable (toList)
import qualified Data.List as List
import Data.Map (Map)
import qualified Data.Map as Map
import Data.Maybe
  ( catMaybes,
    fromJust,
    fromMaybe,
    listToMaybe,
    mapMaybe,
  )
import Data.Tree (Forest, Tree)
import qualified Data.Tree as Tree
import Data.Vector (Vector, (!))
import qualified Data.Vector as Vec
import GHC.Generics (Generic)
import Genomes as G
import LocalBase

-- | Indication of diferent partition types
data PartitionType = MCISP | RMCISP deriving (Show, Read, Enum, Bounded, Eq)

-- | Position of a subgenome X in the indicated genome Y (G or H). The first value (Idx) is the index of X in Y. The second value (Size) is the size of Y minus the size of X. The third value (Ori) indicate whether X and Y are inverted.
data GenomePosition = G !Idx !Size !Ori | H !Idx !Size !Ori deriving (Show, Eq)

instance Orientable GenomePosition where
  getOri (G _ _ LR) = LR
  getOri (G _ _ RL) = RL
  getOri (H _ _ LR) = LR
  getOri (H _ _ RL) = RL

  invOri (G idx remSize LR) = G (coerce remSize - idx + 2) remSize RL
  invOri (G idx remSize RL) = G (coerce remSize - idx + 2) remSize LR
  invOri (H idx remSize LR) = H (coerce remSize - idx + 2) remSize RL
  invOri (H idx remSize RL) = H (coerce remSize - idx + 2) remSize LR

instance Orientable (Maybe GenomePosition) where
    getOri Nothing = LR
    getOri (Just gp) = getOri gp
    invOri = fmap invOri

-- | Tmin is the minimal set of set T, where
--  T is the set of genomes with different
--  weights in g or h
data Tmin = Tmin !PartitionType !Sign HSSTree deriving (Show, Eq)

-- | SuffixTree from the Hitting Set Algorithm, we use it to construct the set Tmin.
--  Leafs correspondent to suffixes started in G are marked with G and Leafs correspondent
--  to suffixes started in H are marked with H.
--  Each node has SumG with the number of leafs of the subtree marked as G and
--  SumH with the number of leafs of the subtree marked as H.
--  Size with the genome size
--  does not contain the characters marking the end of strings.
data HSSTree = HSRoot
  { gSize :: !Size,
    hSize :: !Size,
    rChildren :: [HSSubTree]
  }
  deriving (Eq)

type HSSubTree = Tree HSLabel

instance Show HSSTree where
  show HSRoot {..} =
    "Genome Size:"
      ++ show gSize
      ++ ", "
      ++ show hSize
      ++ "\n"
      ++ Tree.drawForest (map (fmap show) rChildren)

type SumG = Int

type SumH = Int

data HSLabel = HSLabel
  { hsInfo :: !Info,
    hsPref :: !IdxPair
  }
  deriving (Eq)

data IsInG = InG | InH deriving (Eq, Show)

-- A leaf is available if it is proper (the correspondent infix is not empty)
-- and its first intergenic region is not a breakpoint or after a breakpoint
-- the lastBreakDistance is the number of characters after the last breakpoint
-- of the correspondent suffix
data Info
  = LeafInfo
      { lastBreakDistance :: !Int,
        isInG :: !IsInG,
        isRev :: !Ori,
        available :: !Bool
      }
  | NodeInfo
      { sumG :: !SumG,
        sumH :: !SumH
      }
  deriving (Eq)

instance Show HSLabel where
  show HSLabel {..} =
    "Prefix:" ++ show hsPref
      ++ ( case hsInfo of
             LeafInfo {..} ->
               (if isInG == InG then " - G" else " - H")
                 ++ ", lbd:"
                 ++ show lastBreakDistance
                 ++ (if isRev == LR then "" else " R")
                 ++ (if available then "" else " X")
             NodeInfo {..} -> " - sumG: " ++ show sumG ++ ", sumH:" ++ show sumH
         )

data STree
  = Node [(IdxPair, STree)]
  | Leaf
  deriving (Show)

makeTmin :: PartitionType -> Genome -> Genome -> Tmin
makeTmin ptype g h = Tmin ptype (genomeIsSigned g) $ makeHSSTree sTree
  where
    l1 = interleaveListRepresentation g
    l2 = interleaveListRepresentation h
    rl1 = interleaveListRepresentation (invOri g)
    rl2 = interleaveListRepresentation (invOri h)
    str = case ptype of
      MCISP -> Vec.concat [l1, Vec.singleton "$", l2, Vec.singleton "#"]
      RMCISP -> Vec.concat [l1, Vec.singleton "$", l2, Vec.singleton "#", rl1, Vec.singleton "%", rl2, Vec.singleton "&"]
    -- Markers for end of strings
    markers = case ptype of
      MCISP -> ["$", "#"]
      RMCISP -> ["$", "#", "%", "&"]
    getG t = case hsInfo . Tree.rootLabel $ t of
      LeafInfo {..} -> if isInG == InG then 1 else 0
      NodeInfo {..} -> sumG
    getH t = case hsInfo . Tree.rootLabel $ t of
      LeafInfo {..} -> if isInG == InH then 1 else 0
      NodeInfo {..} -> sumH

    -- Construct the Suffix Tree
    sTree = go [0 .. length str - 1]
      where
        go :: [Int] -> STree
        go [] = Leaf
        go ss =
          Node
            [ (IdxPair str (begin_idx -1) end_idx, go ssr)
              | (_, sufs) <- map (second reverse) . Map.toList . suffixMap str $ ss,
                (begin_idx, end_idx, ssr) <- [findEdge sufs]
            ]

        -- input is a list of suffixes and output the begin and end of the edge infix
        -- (first element and last element + 1) and a list of the subtree suffixes.
        findEdge :: [Int] -> (Int, Int, [Int])
        findEdge [] = error patternError
        findEdge [s] = (s, length str, [])
        findEdge sss@(a_idx : ss)
          | null [c_idx | c_idx <- ss, str ! a_idx /= str ! c_idx] =
            let (_, end_idx, ss') =
                  findEdge ((a_idx + 1) : [c_idx + 1 | c_idx <- filter (/= length str - 1) ss])
             in (a_idx, end_idx, ss')
          | otherwise = (a_idx, a_idx, sss)

    -- Convert the suffix tree to HSSTree
    makeHSSTree Leaf = error patternError
    makeHSSTree (Node edges) =
      HSRoot (genomeSize g) (genomeSize h) subTrees
      where
        subTrees = concat $ mapMaybe makeSubHSSTree edges
        makeSubHSSTree e@(idxPair, t) =
          if
              | not $ validBeginIR (getHead idxPair) -> Nothing
              | ipSize idxPair == 1 -> Just $ go False e
              | otherwise ->
                Just $
                  go False (getHeadPair idxPair, Node [(dropHead idxPair, t)])

    go withHead (idxPair_, st) =
      case st of
        Leaf ->
          pure $
            flip Tree.Node [] $
              if
                  | "$" `elem` p__ ->
                    let p = Vec.takeWhile (/= "$") p__
                     in HSLabel
                          (LeafInfo 0 InG LR (getHead idxPair__ /= "$"))
                          (ipTakePrefix idxPair__ (length p))
                  | "#" `elem` p__ ->
                    let p = Vec.takeWhile (/= "#") p__
                     in HSLabel
                          (LeafInfo 0 InH LR (getHead idxPair__ /= "#"))
                          (ipTakePrefix idxPair__ (length p))
                  | "%" `elem` p__ ->
                    let p = Vec.takeWhile (/= "%") p__
                     in HSLabel
                          (LeafInfo 0 InG RL (getHead idxPair__ /= "%"))
                          (ipTakePrefix idxPair__ (length p))
                  | "&" `elem` p__ ->
                    let p = Vec.takeWhile (/= "&") p__
                     in HSLabel
                          (LeafInfo 0 InH RL (getHead idxPair__ /= "&"))
                          (ipTakePrefix idxPair__ (length p))
                  | otherwise -> error patternError
        (Node edges) ->
          if size == 0
            then subTrees
            else [Tree.Node (HSLabel (NodeInfo sG sH) idxPair) subTrees]
          where
            size = ipSize idxPair
            sG = sum . map getG $ subTrees
            sH = sum . map getH $ subTrees
            subTrees = concatMap makeSubHSSTree edges
            makeSubHSSTree = go withHead'
            (withHead', idxPair) =
              ( \h ->
                  if validEndIR h
                    then (False, idxPair__)
                    else (True, dropLast idxPair__)
              )
                $ getLast idxPair__
      where
        p__ = ipSlice idxPair__
        idxPair__ = if withHead then addHead idxPair_ else idxPair_

class WalkDownInfo w where
  goDown :: w -> HSSubTree -> [(w, HSSubTree)]

-- information stored while walking on the tree
data SearchDownInfo = SearchDown
  { subStr :: [IdxPair], -- node's substring
    accPref :: [IdxPair] -- accummulated prefix
  }

idxsToVector :: [IdxPair] -> Vector ByteString
idxsToVector = Vec.concat . map ipSlice . reverse

makeSearchDownInfo :: HSSubTree -> SearchDownInfo
makeSearchDownInfo t@(Tree.Node HSLabel {..} children) =
  if ipSize hsPref >= 2
    then SearchDown [ipTakePrefix hsPref 2] [hsPref]
    else SearchDown [ipTakePrefix hsPref 1] [hsPref]

instance WalkDownInfo SearchDownInfo where
  goDown SearchDown {..} currentNode =
    let (Tree.Node HSLabel {..} children) = currentNode
     in case hsInfo of
          LeafInfo {..} -> []
          NodeInfo {..} -> fmap aux children
    where
      aux t@(Tree.Node HSLabel {..} children) =
        (,t) $
          SearchDown
            (ipTakePrefix hsPref 2 : accPref)
            (hsPref : accPref)

data UpdateInfo = Update
  { currentStrs :: Map ByteString [Int], -- remaining of selected suffixes
    currentGenome :: Vector ByteString, -- genome with breakpoint
    charDistances :: [(Idx, Idx)], -- number of characters until the breakpoint is reached and number of characters until update starts
    leafInG :: IsInG -- whether the suffix's leaf is in G
  }
  deriving (Show)

makeUpdateInfo :: Map ByteString [Int] -> Vector ByteString -> UpdateInfo
makeUpdateInfo m v = Update m v [] InG

instance WalkDownInfo (Maybe UpdateInfo) where
  goDown Nothing _ = []
  goDown (Just up@Update {..}) currentNode =
    case hsInfo of
      LeafInfo {..} -> []
      NodeInfo {..} -> map aux children
    where
      (Tree.Node HSLabel {..} children) = currentNode
      aux :: HSSubTree -> (Maybe UpdateInfo, HSSubTree)
      aux t@(Tree.Node HSLabel {..} children)
        | ipSize hsPref == 0 = (Just up, t)
        | notFound = (Nothing, t)
        | otherwise = (Just $ up {currentStrs = sufMap'}, t)
        where
          (sufMap', notFound) = moveSuffixMap currentGenome hsPref currentStrs

breakNode :: HSSubTree -> Idx -> IsInG -> HSSubTree
breakNode (Tree.Node HSLabel {..} children) idx inG =
  case hsInfo of
    LeafInfo {..} ->
      Tree.Node
        ( HSLabel
            ( if
                  | not available -> NodeInfo 0 0
                  | isInG == InG -> NodeInfo 1 0
                  | otherwise -> NodeInfo 0 1
            )
            pl
        )
        [Tree.Node (HSLabel hsInfo {available = False} pr) children]
    NodeInfo {..} ->
      Tree.Node
        (HSLabel (NodeInfo sumG sumH) pl)
        [ Tree.Node
            ( if inG == InG
                then HSLabel (NodeInfo (sumG - 1) sumH) pr
                else HSLabel (NodeInfo sumG (sumH - 1)) pr
            )
            children
        ]
  where
    (pl, pr) = ipSplitAt hsPref (coerce idx - 1)

updateNode :: UpdateInfo -> Idx -> HSSubTree -> (UpdateInfo, HSSubTree)
updateNode up@Update {..} bDist t@(Tree.Node hsl children) =
  case hsInfo hsl of
    info@LeafInfo {..} ->
      let lastBreakDistance' = max (coerce bDist) lastBreakDistance
          bd = bDist - (coerce . ipSize $ hsPref hsl)
          ud = coerce lastBreakDistance - (coerce . ipSize $ hsPref hsl)
       in ( up
              { charDistances =
                  [ (bd, ud)
                    | bd > -2
                  ],
                leafInG = isInG
              },
            (if ud <= 0 && bd <= -2 then (\node -> breakNode node (- bd) leafInG) else id) $
              Tree.Node
                ( hsl
                    { hsInfo =
                        info
                          { lastBreakDistance = lastBreakDistance',
                            available = available && bd <= -2
                          }
                    }
                )
                children
          )
    NodeInfo {} -> (up {charDistances = charDistances'}, t')
      where
        ((t', _), charDistances') = second catMaybes $ List.mapAccumL aux (t, 0) charDistances
        aux :: (HSSubTree, Idx) -> (Idx, Idx) -> ((HSSubTree, Idx), Maybe (Idx, Idx))
        aux (t@(Tree.Node hsl@HSLabel {..} children), decVal) (bd_, ud_) =
          case hsInfo of
            LeafInfo {} -> error patternError
            NodeInfo {..} ->
              let bd = bd_ - decVal
                  ud = ud_ - decVal
                  bd' = bd - (coerce . ipSize $ hsPref)
               in if
                      | ud <= 0 && bd' > -2 ->
                        if leafInG == InG
                          then
                            ( ( Tree.Node
                                  hsl
                                    { hsInfo = hsInfo {sumG = sumG - 1}
                                    }
                                  children,
                                decVal
                              ),
                              Just (bd', ud)
                            )
                          else
                            ( ( Tree.Node
                                  hsl
                                    { hsInfo = hsInfo {sumH = sumH - 1}
                                    }
                                  children,
                                decVal
                              ),
                              Just (bd', ud)
                            )
                      | bd > 0 && ud <= 0 && bd' <= -2 ->
                        ((breakNode t (- bd') leafInG, decVal + bd + 1), Nothing)
                      | otherwise -> ((t, decVal), Just (bd', ud - (coerce . ipSize $ hsPref)))

-- | Recover a sub-genome of some block also returns
--  a genome position of one of its occurrences
getGenome :: Tmin -> Maybe (Genome, GenomePosition)
getGenome (Tmin ptype sign HSRoot {..}) =
  fmap toGenome . safeMinimum
    . mapMaybe (\t -> walkDown (makeSearchDownInfo t) t)
    $ rChildren
  where
    toGenome :: ([IdxPair], Size, Size, IsInG, Ori) -> (Genome, GenomePosition)
    -- Size is the size of the suffix after the sub-genome occurrence
    -- and the boolean indicate whether the occurrence is in A
    toGenome (idxPairs, _, suf_size_, inG, ori) =
      (x,) $
        if inG == InG
          then G idx (gSize - genomeSize x) LR
          else H idx (hSize - genomeSize x) LR
      where
        -- Sum 1 to correct for the position in the case where there is only one gene in the genome
        suf_size = if genomeSize x_ > 1 then suf_size_ else suf_size_ + 1
        x_ = interleaveListToGenome (idxsToVector idxPairs) sign
        (idx, x) = case ori of
          LR ->
            ( coerce $
                (if inG == InG then gSize else hSize) - suf_size `div` 2 - genomeSize x + 1,
              x_
            )
          RL -> (coerce $ suf_size `div` 2 + 1, invOri x_)

    -- Find minimum element of T'.
    walkDown :: SearchDownInfo -> HSSubTree -> Maybe ([IdxPair], Size, Size, IsInG, Ori)
    walkDown gd@SearchDown {..} currentNode =
      let (Tree.Node HSLabel {..} children) = currentNode
       in case hsInfo of
            LeafInfo {..} ->
              if available
                then Just (subStr, sum . map ipSize $ subStr, ipSize hsPref - 2, isInG, isRev)
                else Nothing
            NodeInfo {..} ->
              if sumG /= sumH
                then
                  let inG = if sumG > sumH then InG else InH
                   in (\(ori, sp) -> (subStr, sum . map ipSize $ subStr, sp, inG, ori))
                        . second (subtract 2)
                        <$> sum_pfx inG 0 currentNode
                else safeMinimum . mapMaybe (uncurry walkDown) $ goDown gd currentNode
    sum_pfx :: IsInG -> Size -> HSSubTree -> Maybe (Ori, Size)
    sum_pfx inG acc (Tree.Node HSLabel {..} children) =
      case hsInfo of
        LeafInfo {..} ->
          let acc' = acc + ipSize hsPref
           in if inG == isInG && lastBreakDistance < coerce acc' - 1
                then Just (isRev, acc')
                else Nothing
        NodeInfo {..} -> listToMaybe $ mapMaybe (sum_pfx inG (ipSize hsPref + acc)) children
    safeMinimum [] = Nothing
    safeMinimum l = Just $ List.minimumBy (\(_, a, _, _, _) (_, b, _, _, _) -> a `compare` b) l

-- | Aditional step when updating Tmin if the update is the deletion of a gene.
deleteGene :: GenomePosition -> Genome -> Tmin -> Tmin
deleteGene gp x (Tmin ptype sign root@HSRoot {..}) =
  Tmin ptype sign $ root {rChildren = map aux rChildren}
  where
    v = interleaveListRepresentation x
    aux t@(Tree.Node l@HSLabel {..} children) =
      if not (IdxPair v 0 (Vec.length v) `ipIsPrefixOf` hsPref)
        then t
        else case hsInfo of
          info@LeafInfo {..} -> Tree.Node l {hsInfo = hsInfo {available = False}} children
          info@NodeInfo {..} ->
            case gp of
              G {} -> Tree.Node l {hsInfo = hsInfo {sumG = sumG}} children
              H {} -> Tree.Node l {hsInfo = hsInfo {sumH = sumH}} children

-- | updateTmin with a new breakpoint, receives the pair of genomes, a genome position indicating the position of the duo representing the breakpoint (position of a genome with two elements), and the Tmin.
updateTmin :: Genome -> Genome -> (Maybe Genome, GenomePosition) -> Tmin -> Tmin
updateTmin g h (delM, breakPos) (Tmin ptype sign root@HSRoot {..}) =
  Tmin ptype sign $ root {rChildren = map aux rChildren}
  where
    (prefSize, k, inG, ori, n) = case breakPos of
      (G gidx _ ori) -> (gidx, g, InG, ori, genomeSize g)
      (H gidx _ ori) -> (gidx, h, InH, ori, genomeSize h)
    bDist = 2 * (coerce n - prefSize) - 1
    v = interleaveListRepresentation k
    sufMap = suffixMap v . take (coerce prefSize) . evens $ [0 ..]
    vdelM = (\v -> IdxPair v 0 (Vec.length v)) . interleaveListRepresentation <$> delM

    -- search suffixes in each subtree of root
    -- Note: a first node never have a breakpoint
    aux t@(Tree.Node l@HSLabel {..} children) =
      ( case vdelM of
          Just idxPairX -> \t'@(Tree.Node l'@HSLabel {..} children') ->
            -- In the case of a deletion we must decrement one child of the root
            if not (idxPairX `ipIsPrefixOf` hsPref)
              then t'
              else case hsInfo of
                info@LeafInfo {..} -> Tree.Node l' {hsInfo = hsInfo {available = False}} children'
                info@NodeInfo {..} ->
                  case breakPos of
                    G {} -> Tree.Node l' {hsInfo = hsInfo {sumG = sumG - 1}} children'
                    H {} -> Tree.Node l' {hsInfo = hsInfo {sumH = sumH - 1}} children'
          Nothing -> id
      )
        $ if notFound
          then t
          else snd . go (makeUpdateInfo sufMap' v) $ t
      where
        (sufMap', notFound) = moveSuffixMap v hsPref sufMap

    -- return the updated node
    go :: UpdateInfo -> HSSubTree -> (Maybe UpdateInfo, HSSubTree)
    go info0 t0@(Tree.Node l@HSLabel {..} children) =
      case goDown (Just info0) t0 of
        [] -> case hsInfo of
          LeafInfo {..} ->
            let (info', t'@(Tree.Node l' _)) = updateNode info0 bDist t0
             in if inG == isInG && Map.null (currentStrs info0) && ori == isRev
                  then (Just info', t')
                  else (Nothing, t0)
          NodeInfo {..} -> (Nothing, t0)
        infosTrees ->
          let (infos_, children') =
                unzip $
                  map
                    ( \(mInfo, t) ->
                        case mInfo of
                          Nothing -> (Nothing, t)
                          Just info -> go info t
                    )
                    infosTrees
              infos = catMaybes infos_
              charDistances' = concatMap charDistances infos
              info' = (head infos) {charDistances = charDistances'}
              (info'', t') =
                (if null charDistances' then (info',) else updateNode info' bDist) $
                  Tree.Node l children'
           in if null infos
                then (Nothing, t0)
                else (Just info'', t')

data IdxPair = IdxPair (Vector ByteString) !Int !Int

instance Show IdxPair where
  show idxPair =
    " ["
      ++ ( List.intercalate ", "
             . map (\p -> let sp = BS.unpack p in (if tail sp == show (maxBound :: Int) then head sp : "inf" else sp))
             . Vec.toList
             $ ipSlice idxPair
         )
      ++ "]"

-- show idxPair = List.intercalate ", " . map (\p -> let sp = show p in (if tail sp == show (maxBound :: Int) then head sp : "inf" else sp) $ ipSlice idxPair

instance Eq IdxPair where
  (IdxPair _ lidx1 ridx1) == (IdxPair _ lidx2 ridx2) = lidx1 == lidx1 && ridx2 == ridx2

addHead :: IdxPair -> IdxPair
addHead (IdxPair v lidx ridx) = IdxPair v (lidx - 1) ridx

getHead :: IdxPair -> ByteString
getHead (IdxPair v lidx _) = v ! lidx

getLast :: IdxPair -> ByteString
getLast (IdxPair v _ ridx) = v ! (ridx - 1)

dropHead :: IdxPair -> IdxPair
dropHead (IdxPair v lidx ridx) = IdxPair v (lidx + 1) ridx

dropLast :: IdxPair -> IdxPair
dropLast (IdxPair v lidx ridx) = IdxPair v lidx (ridx - 1)

getHeadPair :: IdxPair -> IdxPair
getHeadPair (IdxPair v lidx ridx) = IdxPair v lidx (lidx + 1)

ipSize :: IdxPair -> Size
ipSize (IdxPair v lidx ridx) = Size (ridx - lidx)

ipSlice :: IdxPair -> Vector ByteString
ipSlice (IdxPair v lidx ridx) = Vec.slice lidx (ridx - lidx) v

ipTakePrefix :: IdxPair -> Int -> IdxPair
ipTakePrefix (IdxPair v lidx ridx) prf = IdxPair v lidx (lidx + prf)

ipDropPrefix :: IdxPair -> Int -> IdxPair
ipDropPrefix (IdxPair v lidx ridx) prf = IdxPair v (lidx + prf) ridx

ipCombine :: IdxPair -> IdxPair -> IdxPair
ipCombine (IdxPair v lidx1 ridx1) (IdxPair _ lidx2 ridx2) =
  assert (ridx1 == lidx2) (IdxPair v lidx1 ridx2)

ipIsPrefixOf :: IdxPair -> IdxPair -> Bool
ipIsPrefixOf idxPair1 idxPair2 = Vec.length v1 <= Vec.length v2 && v1 == Vec.unsafeTake (Vec.length v1) v2
  where
    v1 = ipSlice idxPair1
    v2 = ipSlice idxPair2

ipSplitAt :: IdxPair -> Int -> (IdxPair, IdxPair)
ipSplitAt idxPair i = (ipTakePrefix idxPair i, ipDropPrefix idxPair i)

suffixMap :: Vector ByteString -> [Int] -> Map ByteString [Int]
suffixMap v = List.foldl' step Map.empty
  where
    step m suf_idx =
      if suf_idx >= Vec.length v
        then m
        else Map.alter (f (suf_idx + 1)) (v ! suf_idx) m
    f i Nothing = Just [i]
    f i (Just is) = Just (i : is)

-- | Remove a prefix of all suffixes from the map. The boolean indicates whether there is no suffixes with the given prefix.
moveSuffixMap :: Vector ByteString -> IdxPair -> Map ByteString [Int] -> (Map ByteString [Int], Bool)
moveSuffixMap v pref m = (sufMap', null sufs)
  where
    sufMap' = suffixMap v sufs
    sufs =
      [ s + coerce (ipSize pref - 1)
        | s <- Map.findWithDefault [] a m,
          dropHead pref `ipIsPrefixOf` IdxPair v s (Vec.length v)
      ]
    a = getHead pref
