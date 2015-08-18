
module PCA
    ( fitPCA, transformPCA
    ) where

import PCA.QModels
import Numeric.LinearAlgebra.Data
       (Matrix, Vector, toRows, fromRows)
import qualified Numeric.LinearAlgebra.HMatrix as H


-- | Theta parameteres which we will you for predict in model
data PCAModel = PCAModel
    { scale :: Matrix Double  -- ^ W
    , shift :: Vector Double  -- ^ \mu
    }


fitPCA :: Matrix Double -> PCAModel
fitPCA t =
    let (w,mu) = calculateQ t
    in PCAModel
       { scale = w
       , shift = mu
       }


transformPCA :: PCAModel -> Matrix Double -> Matrix Double
transformPCA (PCAModel w mu) t =
    fromRows
        (map
             (\tn ->
                   tn + mu)
             (toRows (t H.<> w)))
