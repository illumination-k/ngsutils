import requests 
from logging import getLogger


logger = logging.getLogger(__name__)


def marchantia_genome_loader(save=False):
    URL = "http://marchantia.info/download/download/JGI_3.1.fasta.gz"
    logger.info("Downloading http://marchantia.info/download/download/JGI_3.1.fasta.gz")
