import json
import math
from datetime import datetime, timedelta

import ndex2
import pandas as pd
from django.utils import timezone
from pybiomart import Dataset

from apps.clustering.models import PpiNetworkCache


# TODO: Rewrite
def read_ndex_file_4(fn):
    """
    Given an input string/file, parse the network and return two-column array with interaction partners
    :param fn: Imput NDEx file as string
    :return: Printed as strings, two-column array with interaction partners
    """
    lines6 = ""
    # read edges and nodes into arrays
    if ("edges" in fn.split("nodes")[1]):
        lines5 = fn.split("{\"nodes\":[")
        lines3 = lines5[1].split("{\"edges\":[")[0]
        # remove "cyTableColumn" from array containing edges
        if ("cyTableColumn" in lines5[1]):
            lines4 = lines5[1].split("{\"edges\":[")[1].split("{\"cyTableColumn\":[")[0]
            lines4 = lines4[:-4]
        # take protein name from networkAttributes or nodeAttributes if it is defined there.
        elif ("networkAttributes" in lines5[1]):
            lines4 = lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[0]
            lines4 = lines4[:-4]
            if ("nodeAttributes" in lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[
                1] and "UniprotName" in lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[1]):
                lines6_temp = \
                    lines5[1].split("{\"edges\":[")[1].split("{\"networkAttributes\":[")[1].split(
                        "{\"nodeAttributes\":[")[
                        1]
                lines6 = lines6_temp.split("{\"edgeAttributes\":[")[0]
        else:
            lines4 = lines5[1].split("{\"edges\":[")[1]
    # check if edge-array comes before node-array in file
    elif ("edges" in fn.split("nodes")[0]):
        lines5 = fn.split("{\"nodes\":[")
        lines3 = lines5[1].split("]},")[0] + "]]]"
        lines4 = lines5[0].split("{\"edges\":[")[1][:-4]
    # lines3 contains the nodes, lines4 the edges, lines6 contains nodeAttributes (information from the ndex file usable for the conversion from node IDs to gene IDs)
    # remove signs to allow automatic json to array conversion
    lines3.replace("@", "")
    lines3.replace("uniprot:", "uniprot")
    lines3.replace("signor:", "signor")
    lines3.replace(" ", "")
    lines3.replace("ncbigene:", "")
    lines3.replace("\\n", "")
    lines33 = lines3[:-3].replace("}]", "")
    node_line = lines33.replace("ncbigene:", "")
    nodelinesplit = node_line.split(", ")
    dictlist = []
    # node dict is later filled with keys (node IDs) and the values are NCBI gene IDs
    node_dict = {}
    if not (node_line.endswith("}")):
        node_line = node_line + "}"
    node_line_2 = "[" + node_line + "]"
    tmp2 = json.loads(node_line_2)
    node_dict_2 = {}
    # iterate over lines in nodeAttributes
    if not (lines6 == ""):
        lines6 = "[" + lines6
        # get array with nodeAttributes for current line
        tmp4 = json.loads(lines6[:-4])
        # if node element has attribute "GeneName_A", then the NCBI ID is given in the nodeAttributes
        for item in tmp4:
            if (item['n'] == "GeneName_A"):
                # use node ID and NCBI ID
                node_dict_2[item['po']] = item['v']
    # print(str(item['po']) + " " + str(item['v']))
    # print(node_dict_2)
    dataset = Dataset(name='hsapiens_gene_ensembl',
                      host='http://www.ensembl.org')
    conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'])
    conv_genelist = conv['Gene name'].tolist()
    for item in tmp2:
        dictlist.append(item)
        # write conversion from node ID to gene ID in dictionary, based on nodeAttributes from the data
        if ('r' in item):
            if (any(c.islower() for c in item['r'])):
                gene_name = item['n']
                if (gene_name in conv_genelist):
                    gene_nbr = conv.index[conv['Gene name'] == gene_name]
                    gene_nbr1 = conv.loc[gene_nbr, 'NCBI gene ID'].values[0]
                    node_dict[item['@id']] = gene_nbr1
            # print(item)
            else:
                node_dict[item['@id']] = item['r']
        # print(item)
        else:
            if (item['n'].isdigit()):
                # if gene ID is in node attributes
                # print(item)
                node_dict[item['@id']] = item['n']
            elif (item['n'] in node_dict_2):
                # otherwise use conversion table to convert gene ID to NCBI ID
                gene_name = node_dict_2[item['n']]
                # print(gene_name)
                if (gene_name in conv_genelist):
                    gene_nbr = conv.index[conv['Gene name'] == gene_name]
                    gene_nbr1 = conv.loc[gene_nbr, 'NCBI gene ID'].values[0]
                    node_dict[item['@id']] = gene_nbr1
        # print(gene_nbr1)
    # print(node_dict)
    # remove signs from string to allow json conversion
    lines4.replace("@", "")
    lines4.replace("uniprot:", "uniprot")
    lines4.replace("signor:", "signor")
    lines4.replace(" ", "")
    lines4 = lines4.replace("]", "")
    edge_line = lines4.rstrip()
    edge_line_2 = "[" + edge_line + "]"
    edgelinesplit = edge_line.split(", ")
    edgelist = []
    tmp4 = json.loads(edge_line_2)
    # get dictionary with gene names and NCBI IDs (entrezgene_id)
    dataset = Dataset(name='hsapiens_gene_ensembl',
                      host='http://www.ensembl.org')
    conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'])
    ret = []
    # convert node IDs in edges to NCBI IDs
    for item in tmp4:
        # print(item)
        if (item['s'] in node_dict and item['t'] in node_dict):
            source = node_dict[item['s']]
            target = node_dict[item['t']]
            # print(source)
            # print(target)
            if (source != target and not (math.isnan(float(source))) and not (math.isnan(float(target)))):
                baz = [str(int(source)), str(int(target))]
                ret.append("\t".join(baz))
    # print("\n".join(ret))
    return ("\n".join(ret))


def import_ndex(network_id, force_update=False):
    """
    Download and process the PPI network directly from ndexbio.org
    :param network_id: String, UUID of the network to download
    :param force_update: Boolean, if true the cached version will be ignored and updated
    :return: String, one line per interaction, seperated by tabs
    """
    ndex_server = 'public.ndexbio.org'

    # --- Check if we can use a cached version
    # Connect to NDEx server anonymously, download metadata and get modification time
    network_metadata = ndex2.client.Ndex2(ndex_server) \
        .get_network_summary(network_id)
    network_modification_time = datetime.fromtimestamp(network_metadata['modificationTime'] / 1000.0, tz=timezone.utc)

    # Try and retrieve a cached version. Check if the modification date is within spec, return the cached network
    if not force_update:
        try:
            ppi_network_cache = PpiNetworkCache.objects.get(network_id=network_id)
            datetime_now = timezone.now()
            # The network data modification date must be the same as the one just retrieved,
            # the network cache must have been created within the last 24h
            if ppi_network_cache.data_last_modified == network_modification_time and \
                    datetime_now - timedelta(hours=24) <= ppi_network_cache.last_modified:
                print(f'Network cached on {ppi_network_cache.last_modified.isoformat()}')
                return ppi_network_cache.network_string
        except PpiNetworkCache.DoesNotExist:
            # Download and generate network if no cache exists
            pass

    # Import NDEx from server based on UUID
    nice_cx_network = ndex2.create_nice_cx_from_server(server=ndex_server, uuid=network_id)

    # --- Create a node_id to gene_id dict which maps from the node_id to the gene_id
    node_to_gene_df = pd.DataFrame([x[1] for x in nice_cx_network.get_nodes()]) \
        .rename({'@id': 'Node ID', 'n': 'Gene name'}, axis='columns')

    # If we are using APID, then we need to use another attribute
    if network_id == '9c38ce6e-c564-11e8-aaa6-0ac135e8bacf':
        node_to_gene_df['Gene name'] = node_to_gene_df['Node ID'].map(
            lambda x: nice_cx_network.get_node_attribute_value(x, 'GeneName_A'))

    # --- Create gene_id to other_id dict which maps from gene_id to other ID e.g NCBI IDs
    query_attributes = ['external_gene_name', 'entrezgene_id']
    gene_mapping_df = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org').query(
        attributes=query_attributes).dropna()
    gene_mapping_df.columns = query_attributes
    # set the Gene name (the one used in the networks as ID). Then convert
    # the entrez IDs into int and then to string
    gene_mapping_df = gene_mapping_df \
        .drop_duplicates(subset=['external_gene_name'], keep='first') \
        .set_index('external_gene_name') \
        .astype(int).astype(str)
    # Create the mapping dict
    gene_mapping_dict = gene_mapping_df.to_dict()['entrezgene_id']  # Get the entrez IDs

    # --- Apply gene mapping from gene name to NCBI IDs to the note_to_gene_df and drop missing values
    node_to_gene_df['Gene name'] = node_to_gene_df['Gene name'].map(gene_mapping_dict)
    node_to_gene_dict = node_to_gene_df \
        .set_index('Node ID') \
        .dropna() \
        .to_dict()['Gene name']

    # --- Create the network PPI file
    # Iterate over all edges
    result_list = []
    for _, edge in nice_cx_network.get_edges():
        edge_source = edge.get('s')
        edge_target = edge.get('t')
        if edge_source != edge_target:
            # Convert source and target to NCBI IDs and write into string
            try:
                result_list.append(node_to_gene_dict[edge_source] + '\t' + node_to_gene_dict[edge_target])
            except KeyError:
                # If no mapping can be found, skip this node
                continue

    # --- Save version to cache (db) and return result network string
    result_string = '\n'.join(result_list)
    PpiNetworkCache.objects.update_or_create(
        network_id=network_id,
        defaults={'data_last_modified': network_modification_time, 'network_string': result_string}
    )
    return result_string
