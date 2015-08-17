import Test.Tasty

import qualified QModels as QM

main :: IO ()
main = defaultMain tests

tests :: TestTree
tests = testGroup "PCA Tests" [QM.tests]
