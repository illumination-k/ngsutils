from logging import getLogger
import anndata 
import pandas as pd


logger = getLogger(__name__)


def read_featureCounts(path_to_featureCounts_data: str, path_to_conditions: str=None) -> anndata.core.anndata.AnnData:
    """ read data from featureCounts output

    Args:
        path_to_featureCounts_data: path to featureCounts file
        path_to_conditions: path to condition file (csv format is supported, index is sample id and first columns are sample condition)
    Returns:
        anndata.core.anndata.AnnData: AnnData object with Chr, Start, End, Strand, Length in var and conditions in obs

    """
    featureCounts_df = pd.read_csv(path_to_featureCounts_data, sep="\t", skiprows=1, index_col=0)
    if path_to_conditions is not None:
        conditions_df = pd.read_csv(path_to_conditions, sep=",", index_col=0)
        adata = anndata.AnnData(featureCounts_df.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1),
                                obs=featureCounts_df[['Chr', 'Start', 'End', 'Strand', 'Length']]).T
        adata.obs = conditions_df
    else:
        print("WARNING: If you want to DEG analysis, please load conditions")
        adata = anndata.AnnData(featureCounts_df.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1),
                                obs=featureCounts_df[['Chr', 'Start', 'End', 'Strand', 'Length']]).T
    return adata