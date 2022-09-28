import GenomesCheck as GC
import PartitionCheck as PC
import BalancerCheck as BC

main :: IO ()
main = do
  ans <- fmap and . sequence $
    [ return True
    , GC.tests
    , PC.tests
    , BC.tests
    ]
  print ans
