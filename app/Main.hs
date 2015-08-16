{-# LANGUAGE GADTs #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TemplateHaskell #-}
module Main where

import Frames
import Frames.CSV (rowGen, RowGen(..), readTableOpt)
import Pipes (Producer)
import qualified Data.Foldable as F


tableTypes'
    rowGen
    { rowTypeName = "Iris"
    , columnNames = [ "sepal length"
                    , "sepal width"
                    , "petal length"
                    , "petal width"
                    , "iris name"]
    }
    "data/iris.csv"

irisStream :: Producer Iris IO ()
irisStream = readTableOpt irisParser "data/iris.csv"

loadIris :: IO (Frame Iris)
loadIris = inCoreAoS irisStream

restructureIris :: ( CanDelete IrisName rs
                   , rs' ~ RDelete IrisName rs
                   , AllAre Double (UnColumn rs')
                   , AsVinyl rs') => Record rs -> Record (IrisName ': rs')
restructureIris r = frameCons (rget' irisName' r) (rdel [pr|IrisName|] r)

splitXY :: (AllAre Double (UnColumn rs), AsVinyl rs)
        => Record (s :-> Text ': rs) -> (Text, [Double])
splitXY (recUncons -> (h, t)) = (h, recToList t)

loadDataSet :: IO ([Text], [[Double]])
loadDataSet =
    loadIris >>=
    \ds ->
         return $ unzip $ F.foldMap ((: []) . splitXY . restructureIris) ds


main :: IO ()
main = undefined
