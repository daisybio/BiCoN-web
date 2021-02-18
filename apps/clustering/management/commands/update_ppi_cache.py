from tqdm import tqdm

from django.core.management.base import BaseCommand

from apps.clustering.tasks.ndex_processing import import_ndex

network_id_dict = {
    'APID': '9c38ce6e-c564-11e8-aaa6-0ac135e8bacf',
    'STRING': '275bd84e-3d18-11e8-a935-0ac135e8bacf',
    'BioGRID': 'becec556-86d4-11e7-a10d-0ac135e8bacf',
    'HPRD': '1093e665-86da-11e7-a10d-0ac135e8bacf'
}


def force_update_ppi_cache():
    for network_name, network_id in tqdm(network_id_dict.items()):
        import_ndex(network_id, force_update=True)


class Command(BaseCommand):
    help = 'Force update the NDEx PPI network cache and save to the database'

    def handle(self, *args, **options):
        force_update_ppi_cache()


if __name__ == '__main__':
    force_update_ppi_cache()
