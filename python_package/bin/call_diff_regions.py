#!/usr/bin/env python
########################################################################
## 02/08/2020
## By Xiang Li,
## lux@gwu.edu
## Peng's Lab
## Version.beta
########################################################################
# Usage 
#python ${EXE_PATH} -b ${INPUT_FILE} -c ${INPUT_NAME} -k ${GENE_LIST_FOLDER}/${GENELISTFILE} -l ${GENELISTFILE: :-4} -r ${RESOLUTION} -f ${FRAGMENTSIZE} -g ${GTFFILE} \
#	-w ${WINDOWSIZE} -n ${NORMALIZATION} -t ${REGIONTYPE} -u ${UP_EXTENSION} -d ${DOWN_EXTENSION} -o ${OUTPUTDIR} -p ${Genic_Partition}
########################################################################

import pandas as pd
import numpy as np
import igraph as ig
from scipy import stats
from optparse import OptionParser
import sys
####################################################################################
## FUNCTIONS
### FUNCTION
def Read_Interaction(_PATH_interaction, _resolution, _col_fore, _col_back):
	PATH_interaction=_PATH_interaction
	col_fore = _col_fore
	col_back  = _col_back
	resolution = _resolution
	
	df_interaction = pd.read_csv(PATH_interaction, sep="\t").fillna(0)
	df_interaction = df_interaction[df_interaction.iloc[:,1]!=df_interaction.iloc[:,2]] ### remove self interaction
	df_interaction.loc[:,'#chr']=df_interaction.iloc[:,0].replace('chr','')
	df_interaction.loc[:,'#chr1']=df_interaction.iloc[:,0]
	df_interaction.loc[:,'x1']=df_interaction.iloc[:,1].astype(int)
	df_interaction.loc[:,'x2']=df_interaction.iloc[:,1].astype(int)+int(resolution)
	df_interaction.loc[:,'chr2']=df_interaction.iloc[:,0]
	df_interaction.loc[:,'y1']=df_interaction.iloc[:,2].astype(int)
	df_interaction.loc[:,'y2']=df_interaction.iloc[:,2].astype(int)+int(resolution)

	df_interaction.loc[:,'log_FC'] = np.log2(df_interaction.loc[:,col_fore].replace(0,0.1) / df_interaction.loc[:,col_back].replace(0,0.1) )
	df_interaction = df_interaction.loc[:,['#chr1','x1','x2','chr2','y1','y2','log_FC', col_fore, col_back]]
	return df_interaction

def Convert_Loops_to_Graph(_df_hic, _weight_col):
	## Assign a list of weight ot graph
	
	## loop format: ['#chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'GeneID', 'weight_cols']
	df_bins = Loops_Return_two_bins_no_dup(_df_hic)
	
	## eliminate float in chr
	df_bins['name'] = df_bins['#chr1'].astype(str).str.split(".",expand=True)[0]+':'+df_bins['x1'].astype(int).astype(str)+'-'+df_bins['x2'].astype(int).astype(str)
	Num_vs = len(df_bins.index)
	## Initiation a graph from loops file 
	graph_tem = ig.Graph()
	graph_tem.add_vertices(Num_vs)
	graph_tem.vs["name"] = df_bins.loc[:,'name']
	df_edge = _df_hic.merge(df_bins, on=['#chr1', 'x1', 'x2']).merge(
		df_bins, left_on=['chr2', 'y1', 'y2'], right_on=['#chr1', 'x1', 'x2'])
	graph_tem.add_edges(df_edge.loc[:, ['index_x','index_y']].values)

	for weight in _weight_col:
		if (weight in _df_hic.columns):
			graph_tem.es[weight] = df_edge.loc[:,weight].values
	return graph_tem

def Loops_Return_two_bins_no_dup(df_hic):
	## Associated by promoter
	second_bin_columns = [3,4,5,0,1,2]+list(range(6,len(df_hic.columns),1))
	df_hic=df_hic.append(pd.DataFrame(df_hic.iloc[:, second_bin_columns].values, columns=df_hic.columns),sort=False).sort_index()
	return df_hic.iloc[:,0:3].drop_duplicates().reset_index().drop('index',axis=1).reset_index()

def convert_cluster2bed(df_cluster, usecol):
	df_tem = df_cluster[usecol].str.split(r"\:|-",expand=True)
	df_tem = pd.concat( [df_tem, df_cluster], axis=1)
	if (df_tem.iloc[0,0].find('chr') == -1):
		df_tem[0] = 'chr'+df_tem[0]
	return df_tem

def convert_bin2bed(df_cluster, col_name):
	df_tem = df_cluster[col_name].str.split(r"\:|-",expand=True)
	df_tem = pd.concat( [df_tem, df_cluster], axis=1)
	if (df_tem.iloc[0,0].find('chr') == -1):
		df_tem[0] = 'chr'+df_tem[0]
	return df_tem

def convert_vs2bed(input_graph, col_name):
	## output first 3 columns is standard bed format
	df_tem = pd.DataFrame(data={col_name:input_graph.vs[col_name]})
	df_tem = pd.concat( [df_tem[col_name].str.split(r"\:|-",expand=True),df_tem], axis=1)
	if (df_tem.iloc[0,0].find('chr') == -1):
		df_tem[0] = 'chr'+df_tem[0]
	return df_tem

def convert_graph_vs_to_df(_input_graph):
	df_vs = pd.DataFrame(data= {"degree":_input_graph.degree()})
	for col in _input_graph.vs.attributes():
		df_vs[col] = _input_graph.vs[col]

	return df_vs

def graph_community_multilevel_Blondel(input_graph, cutoff):
	## input graph should have at least one attribute: name
	df_vs = convert_graph_vs_to_df(input_graph)
	_col_vs_name='name'
	if (input_graph.is_weighted()):
		print ("Weighted Graph Cluster")
		structure = input_graph.community_multilevel(weights=input_graph.es['weight'] ,return_levels=False)
	else:
		structure = input_graph.community_multilevel(return_levels=False)
	df_vs['membership'] = structure.membership
	df_vs_cluster_group = df_vs.groupby('membership')
	
	## Rank each cluster by number of bins
	cluster_name=[]
	cluster_num_vertices=[]
	for df_vs_cluster in df_vs_cluster_group:
		df_vs_inside_cluster = Cluster_Filter_by_Denisty(df_vs_cluster[1], _col_vs_name, 'degree', cutoff)
		#df_vs_inside_cluster =df_vs_cluster[1]
		df_cluster_coordiante = df_vs_inside_cluster[_col_vs_name].str.split(r"\:|-",expand=True)
		cluster_coordinate = 'chr'+df_cluster_coordiante.iloc[0,0]+':'+str(df_cluster_coordiante.iloc[:,1].astype(int).min())+'-'+str(df_cluster_coordiante.iloc[:,2].astype(int).max())
		cluster_name.append(cluster_coordinate) ##0: cluster name
		cluster_num_vertices.append(len(df_vs_inside_cluster)) # 1: num_vertices
	
	df_cluster_output = pd.DataFrame(data={'hub_name':cluster_name,'Num_vertices':cluster_num_vertices}).sort_values('Num_vertices', ascending=False)
	return df_cluster_output, df_vs_cluster_group

def Graph_Pagerank(_input_graph):
	input_graph = _input_graph
	input_graph.vs['pagerank'] = input_graph.pagerank(weights=input_graph.es['weight'])
	return input_graph

### allow a gap size of one window

def Stich_Region_Above_Mean(_graph, _resolution, _gap_size):
	resolution=_resolution
	graph_pagerank = Graph_Pagerank(_graph)
	df_vs_graph = convert_graph_vs_to_df(graph_pagerank)
	df_nodes = convert_cluster2bed(df_vs_graph, 'name')
	df_nodes[1] = df_nodes[1].astype(int)
	df_nodes = df_nodes.sort_values(by=1)
	df_nodes = df_nodes[df_nodes['pagerank'] > df_nodes['pagerank'].mean()] ## Only use nodes > mean
	
	## report stich regions
	Report_list=[]
	reg_chr = str(df_nodes.iloc[0,0])
	reg_start= int(df_nodes.iloc[0,1])
	reg_end = int(reg_start)
	
	for bin1 in df_nodes.iloc[:,1].astype(int):
		if (bin1-reg_end)<=_gap_size*resolution:
			reg_end = bin1
		else:
			Report_list.append([reg_chr+':'+str(reg_start)+'-'+str(reg_end+resolution), _gap_size])
			reg_start = bin1
			reg_end = bin1
	Report_list.append([reg_chr+':'+str(reg_start)+'-'+str(reg_end+resolution), _gap_size])
	
	return pd.DataFrame(data=Report_list, columns=['hub_name', 'merge_level'])

def Stich_Region_Above_global_Mean(_graph, _resolution, _gap_size, _mean):
	resolution=_resolution
	df_vs_graph = convert_graph_vs_to_df(_graph)
	df_nodes = convert_cluster2bed(df_vs_graph, 'name')
	df_nodes[1] = df_nodes[1].astype(int)
	df_nodes = df_nodes.sort_values(by=1)
	df_nodes = df_nodes[df_nodes['pagerank'] > _mean] ## Only use nodes > mean
	Report_list=[]
	if (len(df_nodes)>0):
		## report stich regions
		
		reg_chr = str(df_nodes.iloc[0,0])
		reg_start= int(df_nodes.iloc[0,1])
		reg_end = int(reg_start)

		for bin1 in df_nodes.iloc[:,1].astype(int):
			if (bin1-reg_end)<=_gap_size*resolution:
				reg_end = bin1
			else:
				Report_list.append([reg_chr+':'+str(reg_start)+'-'+str(reg_end+resolution), _gap_size])
				reg_start = bin1
				reg_end = bin1
		Report_list.append([reg_chr+':'+str(reg_start)+'-'+str(reg_end+resolution), _gap_size])
	return pd.DataFrame(data=Report_list, columns=['hub_name', 'merge_level'])


def Return_Sorted_Adjacency_Matrix(_graph, _attr):
	
	## Sort by coordinate
	graph_tem = _graph
	attr      =_attr
	idx_name = [int(str(x).split(":")[1].split("-")[0]) for x in graph_tem.vs['name']]
	
	matrix_tem = pd.DataFrame(data=graph_tem.get_adjacency(attribute=attr), columns=idx_name, index=idx_name)
	df_reindex = pd.DataFrame(data={ 'rank': (stats.rankdata(matrix_tem.columns)-1).astype(int)})
	idx_rank = df_reindex.sort_values(by='rank').index
	## reference https://wil.yegelwel.com/cluster-correlation-matrix/
	return matrix_tem.iloc[idx_rank, :].T.iloc[idx_rank, :]

def Pvalue_Rank_Test_Matrix(_matirx):
    matrix_for_test = _matirx
    
    data_test = matrix_for_test.fillna(0).values.flatten() ## flatten 2d into 1D
    if (len(data_test)>10):
        w, pvalue =stats.wilcoxon(data_test, zero_method='zsplit', alternative='greater', correction=True, mode='approx')  
        # “zsplit”: Includes zero-differences in the ranking process and split the zero rank between positive and negative ones.
    else:
        pvalue=1.0
    return round(pvalue, 2)

def Return_Pvalue_For_Given_Graph(_df_region, _resolution, _matrix):
    df_region = _df_region
    df_regionh_bed = convert_cluster2bed(df_region, 'hub_name').sort_values(by=1)
    resolution = _resolution
    matrix_for_test = _matrix
    
    ## convert each region into bins
    idx_regs = []
    for name_stitch in df_region.hub_name:
        region_loc= name_stitch.split(":")[1].split("-")
        idx_reg = []
        for idx in matrix_for_test.index:
            if ((idx>=int(region_loc[0]))&(idx<=int(region_loc[1]))):
                idx_reg.append(idx)
        idx_regs.append(idx_reg)

    pvalue_region= []
    for i in range(len(idx_regs)):
        for j in range(i+1):
            part_matrix_for_test = matrix_for_test.loc[idx_regs[i],:].T.loc[idx_regs[j], :]
            pvalue_tem = Pvalue_Rank_Test_Matrix(part_matrix_for_test)
            pvalue_region.append([df_region.hub_name[i],df_region.hub_name[j],-np.log10(pvalue_tem)])


    return pd.DataFrame(data=pvalue_region, columns=['reg1', 'reg2', '-log10(pvalue)']).sort_values('-log10(pvalue)', ascending=False)


def Main_For_Diff_Regions(_df_hic, _col_fore, _col_back,  _resolution, _gapsize):
	#Create a weight basing on logFC (logFC < 0)
	_gapsize=2
	_df_hic.loc[:,_col_back+'_weight'] = _df_hic.loc[:,_col_back]*_df_hic.log_FC.apply(lambda x: 0 if x > 0 else(1))
	_df_hic.loc[:,_col_fore+'_weight'] = _df_hic.loc[:,_col_fore]*_df_hic.log_FC.apply(lambda x: 1 if x > 0 else(0))
	weight_list= [_col_back+'_weight', _col_fore+'_weight', _col_back, _col_fore]
	input_graph = Convert_Loops_to_Graph(_df_hic, weight_list)
	#######################################################
	input_graph.es['weight'] = input_graph.es[_col_fore+'_weight']
	structure = input_graph.community_multilevel(weights=input_graph.es['weight'], return_levels=True)
	input_graph = Graph_Pagerank(input_graph)
	global_mean = np.mean(input_graph.vs['pagerank'])
	#######################################################

	### Stich according to pagerank locally
	df_out = pd.DataFrame(columns=['reg1', 'reg2', '-log10(pvalue)'])
	i=0
	for graph_tem in structure[0].subgraphs(): # subgraphs() is much faster than using subgraph(idx)
		if (len(graph_tem.vs)>_gapsize+1): ## #_nodes shoud be > _gapsize
			i+=1
			#df_hubs = Stich_Region_Above_Mean(graph_tem, _resolution, _gapsize)  ## approximaltly 0.1s for each graph
			df_hubs = Stich_Region_Above_global_Mean(graph_tem, _resolution, _gapsize, global_mean)  ## approximaltly 0.1s for each graph
			if(len(df_hubs)>0):
				Back_matrix = Return_Sorted_Adjacency_Matrix(graph_tem, _col_back)
				Fore_matrix = Return_Sorted_Adjacency_Matrix(graph_tem, _col_fore)
				Diff_matrix = Fore_matrix-Back_matrix
				df_out = df_out.append(Return_Pvalue_For_Given_Graph(df_hubs, _resolution, Diff_matrix))
	df_out = df_out.sort_values(by='-log10(pvalue)', ascending=False)
	df_out = df_out[df_out['-log10(pvalue)']>5]
	df_out.to_csv(str(len(df_out))+'_'+_col_fore+'_specific_regions.bed', sep='\t', index=None)
	
	
	#######################################################
	input_graph.es['weight'] = input_graph.es[_col_back+'_weight']
	structure = input_graph.community_multilevel(weights=input_graph.es['weight'], return_levels=True)
	input_graph = Graph_Pagerank(input_graph)
	global_mean = np.mean(input_graph.vs['pagerank'])
	#######################################################

	### Stich according to pagerank locally
	df_out = pd.DataFrame(columns=['reg1', 'reg2', '-log10(pvalue)'])
	i=0
	for graph_tem in structure[0].subgraphs(): # subgraphs() is much faster than using subgraph(idx)
		if (len(graph_tem.vs)>_gapsize+1): ## #_nodes shoud be > _gapsize
			i+=1
			#df_hubs = Stich_Region_Above_Mean(graph_tem, _resolution, _gapsize)  ## approximaltly 0.1s for each graph
			df_hubs = Stich_Region_Above_global_Mean(graph_tem, _resolution, _gapsize, global_mean)  ## approximaltly 0.1s for each graph
			if(len(df_hubs)>0):
				Back_matrix = Return_Sorted_Adjacency_Matrix(graph_tem, _col_back)
				Fore_matrix = Return_Sorted_Adjacency_Matrix(graph_tem, _col_fore)
				Diff_matrix = Back_matrix-Fore_matrix
				df_out = df_out.append(Return_Pvalue_For_Given_Graph(df_hubs, _resolution, Diff_matrix))
	df_out = df_out.sort_values(by='-log10(pvalue)', ascending=False)
	df_out = df_out[df_out['-log10(pvalue)']>2]
	df_out.to_csv(str(len(df_out))+'_'+_col_back+'_specific_regions.bed', sep='\t', index=None)
	return None
### End of Visulization
####################################################################################
### FUNCTION
### FUNCTIONS
def main(argv):
	desc="Collect HiC Interaction in txt format, rank interaction change Hub. Input Format should be: #chr	bin1	bin2	Cond1	Cond2"
	parser = OptionParser(description=desc)
	parser.add_option("-i", "--in", action="store", type="string",
			dest="input_path", help="Path to Input HiC file in txt format", metavar="<file>")
	parser.add_option("-f", "--foreground_name", action="store", type="string",
			dest="fore_name", help="Name of condition as foreground.", metavar="<str>")
	parser.add_option("-b", "--background_name", action="store", type="string",
			dest="back_name", help="Name of condition as background.", metavar="<str>")
	parser.add_option("-r", "--resolution", action="store", type="int",
		dest="res", help="Resolution of HiC txt", metavar="<int>")
	parser.add_option("-g", "--gap_size", action="store", type="float",
		dest="gap_size", help="Gap Size to stitch nodes.", metavar="<float>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 5:
		parser.print_help()
		sys.exit(1)
	
	print (" ")
	print ("Here is the Summary of your input.")
	print ("Input Path of HiC file in txt format: %s" % opt.input_path)
	print ("Foreground Condition: %s" % opt.fore_name)
	print ("Background Condition: %s" % opt.back_name)
	print ("Resolution %i" % opt.res)
	print ("Gap Size to stitch nodes: %g" % opt.gap_size)
	print ("End of Summary.")
	print (" ")
	
	## parameters
	PATH_INPUT=opt.input_path
	col_fore = opt.fore_name
	col_back  = opt.back_name
	resolution = opt.res

	_gapsize=opt.gap_size

#### Main 
	df_hic = Read_Interaction(PATH_INPUT, resolution, col_fore, col_back)
	df_out = Main_For_Diff_Regions(df_hic, col_fore, col_back, resolution, _gapsize)

	print(" ")
#### First GeneBoydy



if __name__ == "__main__":
	main(sys.argv)
