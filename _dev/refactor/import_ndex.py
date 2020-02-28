# from the web
import ndex2
import math
import pandas as pd
import numpy as np
from pybiomart import Dataset


def import_ndex(name):
    # import NDEx from server based on UUID (network contains lists with nodes and edges)
    nice_cx_network = ndex2.create_nice_cx_from_server(server='public.ndexbio.org', uuid=name)
    print(nice_cx_network.metadata)

    # node_dict is for the conversion from Node ID to NCBI ID
    node_dict = {}
    dataset = Dataset(name='hsapiens_gene_ensembl',
                      host='http://www.ensembl.org')
    # get list of genes for later conversion
    conv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'])
    conv_genelist = conv['Gene name'].tolist()
    # iterate over all nodes in network
    for node_id, node in nice_cx_network.get_nodes():
        # node has ID and "name"
        current_node = {'id': node_id, 'n': node.get('n')}
        # gene names are stored in the GeneName variable in this network
        if name == "9c38ce6e-c564-11e8-aaa6-0ac135e8bacf":
            # get GeneName for node
            curr_gene_name = nice_cx_network.get_node_attribute_value(node_id, 'GeneName_A')
            if curr_gene_name in conv_genelist:
                # get index in gene conversion list, convert Gene Name to NCBI ID
                gene_nbr = conv.index[conv['Gene name'] == curr_gene_name]
                gene_nbr1 = conv.loc[gene_nbr, 'NCBI gene ID'].values[0]
                # check if NCBI ID was found
                if not (math.isnan(float(gene_nbr1))):
                    node_dict[node_id] = str(int(gene_nbr1))

        else:
            # if gene name is stored in node name
            curr_gene_name = current_node['n']
            if curr_gene_name in conv_genelist:
                gene_nbr = conv.index[conv['Gene name'] == curr_gene_name]
                gene_nbr1 = conv.loc[gene_nbr, 'NCBI gene ID'].values[0]
                if not (math.isnan(float(gene_nbr1))):
                    node_dict[node_id] = str(int(gene_nbr1))

    edgelist = []
    ret = ""
    # iterate over edges
    for edge_id, edge in nice_cx_network.get_edges():
        source = edge.get('s')
        target = edge.get('t')
        if source in node_dict and target in node_dict and source != target:
            # convert source and target to NCBI IDs and write into string
            curr_edge_str = str(node_dict[source]) + "\t" + str(node_dict[target]) + "\n"
            edgelist.append([node_dict[source], node_dict[target]])
            ret = ret + curr_edge_str
    # return tab separated string
    return ret


def import_ndex2(name):
    # Check when the network was last modified and use local cached version if nothing has changed
    # using the ndex client

    # Import NDEx from server based on UUID
    nice_cx_network = ndex2.create_nice_cx_from_server(server='public.ndexbio.org', uuid=name)

    # --- Create a node_id to gene_id dict which maps from the node_id to the gene_id
    node_to_gene_df = pd.DataFrame([x[1] for x in nice_cx_network.get_nodes()], dtype={'int', 'str'})
    node_to_gene_df.columns = ['Node ID', 'Gene name']

    # If we are using Apid, then we need to use another attribute
    if name == '9c38ce6e-c564-11e8-aaa6-0ac135e8bacf':
        node_to_gene_df['Gene name'] = node_to_gene_df['Node ID'].map(lambda x: nice_cx_network.get_node_attribute_value(x, 'GeneName_A'))


    # --- Create gene_id to other_id dict which maps from gene_id to other ID e.g Entrez
    mapping_df = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org').query(
        attributes=['external_gene_name', 'entrezgene_id']).dropna()
    # set the Gene name (the one used in the networks as ID). Then convert
    # The entrez IDs into int and then to string
    mapping_df = mapping_df.set_index('Gene name').astype(int).astype(str)
    # Create the mapping dict
    mapping_dict = mapping_df.to_dict()['NCBI gene ID']  # Get the entrez IDs

    # --- Create the network PPI file
    # Iterate over all edges
    result_list = []
    for edge_id, edge in nice_cx_network.get_edges():
        edge_source = edge.get('s')
        edge_target = edge.get('t')
        if edge_source in mapping_dict and edge_target in mapping_dict and edge_source != edge_target:
            # convert source and target to NCBI IDs and write into string
            result_list.append(mapping_dict[edge_source] + '\t' + mapping_dict[edge_target])

    # return tab separated string with linebreaks
    return '\n'.join(result_list)


if __name__ == '__main__':
    # APID
    result = import_ndex2('9c38ce6e-c564-11e8-aaa6-0ac135e8bacf')
    result = import_ndex('9c38ce6e-c564-11e8-aaa6-0ac135e8bacf')

    with open("outfiles/new.tsv", "w") as file:
        file.writelines(result)
