import pandas as pd
import numpy as np
import re

def make_group_analytes_unique( grouping_file , delimiter='\t' ) :
    uniqe_grouping_file = '/'.join(grouping_file.split('/')[:-1]) + '/unique_'+grouping_file.split('/')[-1]
    with open( uniqe_grouping_file , 'w' ) as of:
        with open( grouping_file ) as input :
            for line in input :
                vline = line.replace('\n','').split(delimiter)
                gid, gdesc, analytes_ = vline[0], vline[1], list(set(vline[2:]))
                nvec = [gid,gdesc] ; [ nvec.append(a) for a in analytes_ ]
                print ( delimiter.join(nvec) , file = of )

def find_compartment( istr ) :
    return ( re.findall( r'\[(.*?)\]', istr ) )

def read_gene_ensemble_conversion(file_name):
    gene2ens = {} ; non_unique = []
    with open(file_name,'r') as infile:
        if sys.version_info[0] < 3:
            infile.next()
        else:
            next(infile)
        for line in infile:
            words = line.strip().split('\t')
            if len(words)==2 :
                ens, gene = words
                if gene in gene2ens:
                    non_unique.append((gene,ens,gene2ens[gene]))
                else :
                    gene2ens[gene] = ens
        if len(non_unique)>0:
            print(' WARNING ' )
            print( 'FOUND ', len(non_unique), ' NON UNIQUE ENTRIES' )
    return gene2ens

def read_conversions(file_name):
    gene2ens = {} ; non_unique = []
    with open(file_name,'r') as infile:
        if sys.version_info[0] < 3:
            infile.next()
        else:
            next(infile)
        for line in infile:
            words = line.strip().split('\t')
            if len(words)==2 :
                ens, gene = words
                if gene in gene2ens:
                    gene2ens[gene].append(ens)
                else :
                    gene2ens[gene] = [ens]
    return gene2ens

def create_synonyms( convert_file , unique_mapping=False ) :
    # CREATE SYNONYMS
    ens2sym,sym2ens = {},{}
    if unique_mapping:
        sym2ens = read_gene_ensemble_conversion( convert_file )
        ens2sym = { v:k for k,v in sym2ens.items() }
    else :
        sym2ens_list = read_conversions( convert_file )
        ens2sym_list = {}
        for s,L in sym2ens_list.items() :
            for e in L:
                if e in ens2sym_list:
                    ens2sym_list.append(s)
                else:
                    ens2sym_list[e]=[s]
        ens2sym = ens2sym_list
        sym2ens = sym2ens_list
    return ( ens2sym , sym2ens )

## PATHWAY SECTION
class Pathway( object ) :
    def get_next(self):
        pass

class GenericPathway( Pathway ) :
    def __init__(   self , path , gene_name_start = "ENSG0" , gene_mapping = None, # GENE ID MAPPING
                    is_own_pathway = False, list_like = False , add_pathway_prefix='',
                    gene_pos=0 , pathway_pos=1 , pathway_desc_pos=3, seperator='\t' ) :
        self.file = path
        self.prefered_genes = gene_name_start
        self.pathways , self.pathway_names = {},{}
        self.sep = seperator
        self.is_own_pathway = is_own_pathway
        self.gene_mapping = gene_mapping
        self.gene_pos = gene_pos
        self.pathway_pos = pathway_pos
        self.pathway_desc_pos = pathway_desc_pos
        self.list_like = list_like
        self.add_pathway_prefix = add_pathway_prefix
        self.replace_pair = None
        self.internal_table = None
        self.generator = None
        self.active_read = False
        self.restrict_id_to = None
        if 'str' in str(type(self.file)):
            self.read_and_parse()
        else :
            self.parse_df()

    def add_pathway_synonyms ( self , synonym_dict , prefix='' ) :
        for pathway,pathway_name,genes in self.get_generator() :
            synonyms = []
            for g in genes :
                if g in synonym_dict or prefix+g in synonym_dict :
                    k = prefix * ( not g in synonym_dict ) + g
                    sdg = synonym_dict[k]
                    if 'list' in str(type(sdg)) :
                        [ synonyms.append(s) for s in sdg ]
                    if 'str' in str(type(sdg)) :
                        synonyms.append(sdg)
            [ self.pathways[pathway].append(s) for s in synonyms ]

    def make_gmt_pathway_file ( self , filename , verbose=False , delimiter='\t', gene_mapping=None ):
        #if 'None' in str(type( self.generator )) :
        #    self.generator = self.get_generator()
        if not gene_mapping is None:
            self.gene_mapping = gene_mapping
        if 'str' in str(type(filename)) :
            if 'str' in str( type( filename ) ) :
                with open( filename, 'w' ) as o_file :
                    for pathway, pathway_name, genes in self.get_generator() :
                        row = list() ; row . append ( pathway )
                        row . append ( pathway_name )
                        genes_loc = genes
                        if 'dict' in str(type(self.gene_mapping)) :
                            genes_loc = [ self.gene_mapping[gene] if gene in self.gene_mapping.keys() else gene for gene in genes ]
                        [ row.append ( gene ) for gene in list(set(genes_loc)) ]
                        row_str = delimiter.join(row)
                        print(row_str)
                        o_file . write( row_str+'\n' )
        else :
            print ( 'MUST SUPPLY A FILENAME' )

    def make_internal_table(self, verbose=False, delimiter='\t', output_file=None ) :
        self.internal_table = list()
        generator = self.get_generator( )
        for pathway, pathway_name, genes in generator :
            row = list() ; row . append ( pathway )
            row .append ( pathway_name )
            [ row .append ( gene ) for gene in genes ]
            row_str = delimiter.join(row)
            self.internal_table.append( row )

    def parse_df(self):
        if not self.is_own_pathway:
            print('ERROR: OPERATION NOT SUPPORTED')
            exit(1)
        sfcv = set( self.file.columns.values )
        if 'gene' in sfcv and 'pathway' in sfcv and 'description' in sfcv:
            print ( self.file.columns )
        else:
            genes = self.file.index.values
        for gene in genes :
            pathway = self.add_pathway_prefix + gene
            self.pathways[pathway] = [gene]
            if not 'None' in str(type(self.gene_mapping)) and gene in self.gene_mapping.keys() :
                self.pathway_names[pathway] = self.gene_mapping[gene]
            else :
                self.pathway_names[pathway] = gene

    def get_generator_from_df(self):
        self.parse_df()
        for key in self.pathways.keys():
            yield(key,self.pathway_names[key],self.pathways[key])

    def read_and_parse(self):
        with open(self.file) as input:
            pathway, pathway_name, genes = "", "", []
            for line in input:
                lspl = line.split('\t')
                if not 'None' in str(type(self.replace_pair)):
                    pathway = ( self.add_pathway_prefix + lspl[self.pathway_pos] ).replace( self.replace_pair[0],self.replace_pair[1] )
                else:
                    pathway = ( self.add_pathway_prefix + lspl[self.pathway_pos] )
                pathway_name = lspl[self.pathway_desc_pos]
                if self.list_like :
                    # LIST LIKE PATHWAY INVENTORY CANNOT HAVE SELF MAPPING
                    # HERE WE ASSUME ALL GENES FOLLOW THE FIRST GENE_POS
                    genes = [ lspl[ir].replace('\n','') for ir in range(self.gene_pos,len(lspl)) ]
                    if not 'None' in str(type(self.gene_mapping)) :
                        renamed_genes = [ self.gene_mapping[gene] if gene in self.gene_mapping.keys() else gene for gene in genes  ]
                        genes = renamed_genes
                    if pathway in self.pathways :
                        [ self.pathways[pathway].append(gene) for gene in genes ]
                    else :
                        self.pathways[pathway] = genes
                        self.pathway_names[pathway] = pathway_name
                else :
                    if not line.startswith(self.prefered_genes) or len(line)==0:
                        continue
                    gene = lspl[ self.gene_pos ]
                    if self.is_own_pathway :
                        if gene in self.pathways :
                            continue;
                        else:
                            pway=self.add_pathway_prefix + gene
                            self.pathways[pway] = [gene]
                            if not 'None' in str(type(self.gene_mapping)) and gene in self.gene_mapping.keys():
                                self.pathway_names[pway] = self.gene_mapping[gene]
                            else:
                                self.pathway_names[pway] = gene
                    else :
                        if not 'None' in str(type(self.gene_mapping)) and gene in self.gene_mapping.keys():
                            gene = self.gene_mapping[gene]
                        if pathway in self.pathways:
                            self.pathways[pathway].append(gene)
                        else:
                            self.pathways[pathway] = [gene]
                            self.pathway_names[pathway] = pathway_name

    def dump_pathways(self):
        return ( self.pathways,self.pathway_names )

    def get_generator( self, verbose=False ):
        if self.active_read :
            if not self.file is None :
                if 'str' in str(type(self.file)):
                    self.read_and_parse()
                else:
                    self.parse_df()
        if verbose :
            print( self.dump_pathways() )
        for key in self.pathways.keys():
            yield( key , self.pathway_names[key] , self.pathways[key] )

class Reactome( GenericPathway ) :
    def __init__(   self , path , gene_name_start = 'ENSG0' ,pathway_desc_pos=3,
                    gene_mapping = None, lexical_restrictions = None, # GENE ID MAPPING
                    is_own_pathway = False, restrict_id_to = None ) :
        self.file = path
        self.prefered_genes = gene_name_start
        self.pathways , self.pathway_names ,self.pathway_compartments = {},{},{}
        self.is_own_pathway = is_own_pathway
        self.gene_mapping = gene_mapping
        self.gene_pos = 0
        self.pathway_pos = 1
        self.pathway_desc_pos = pathway_desc_pos
        self.list_like = False
        self.add_pathway_prefix = ''
        self.replace_pair = None
        self.lexical_restrictions = lexical_restrictions
        self.lexical_restriction_category_pos=None
        self.skipstr = None
        self.pathway_tag = None
        self.restrict_id_to = restrict_id_to
        self.active_read = False
        if 'str' in str(type(self.file)):
            self .read_and_parse(keep_str='ENSG0')
        else:
            self .parse_df()

    def read_and_parse( self , keep_str=None ) :
        with open( self.file ) as input :
            pathway, pathway_name, genes = "", "", []
            for line in input :
                if not self.skipstr is None :
                    if line.startswith(self.skipstr):
                        continue
                if not keep_str is None :
                    if not line.startswith(keep_str):
                        continue
                lspl = line.split('\t')
                if ('[' in line) and (']' in line) :
                    compartment = find_compartment( line )[0]
                else :
                    compartment = ''
                if not 'None' in str(type(self.replace_pair)) :
                    pathway = ( self.add_pathway_prefix + lspl[self.pathway_pos] ).replace( self.replace_pair[0],self.replace_pair[1] )
                else :
                    pathway = ( self.add_pathway_prefix + lspl[self.pathway_pos] )

                if not self.restrict_id_to is None :
                    if not pathway in self.restrict_id_to:
                        continue
                pathway_name = lspl[ self.pathway_desc_pos ]
                if not self.lexical_restrictions is None :
                    if not np.sum( [ int(lr in pathway_name) for lr in self.lexical_restrictions] )>0 : #or len(compartment)<4
                        continue
                if self.lexical_restriction_category_pos is None :
                    lex_restrict = pathway_name
                else :
                    lex_restrict = lspl[ self.lexical_restriction_category_pos ]

                if self.list_like :
                    # LIST LIKE PATHWAY INVENTORY CANNOT HAVE SELF MAPPING
                    # HERE WE ASSUME ALL GENES FOLLOW THE FIRST GENE_POS
                    genes = [ lspl[ir].replace('\n','') for ir in range(self.gene_pos,len(lspl)) ]
                    if not 'None' in str(type(self.gene_mapping)) :
                        renamed_genes = [ self.gene_mapping[gene] if gene in self.gene_mapping.keys() else gene for gene in genes  ]
                        genes = renamed_genes
                    if pathway in self.pathways :
                        [ self.pathways[pathway].append(gene) for gene in genes ]
                    else :
                        self.pathways[pathway] = genes
                        self.pathway_names[pathway] = pathway_name
                        self.pathway_compartments[pathway] = compartment
                else :
                    if not line.startswith(self.prefered_genes) or len(line)==0:
                        continue
                    gene = lspl[ self.gene_pos ]
                    if self.is_own_pathway :
                        if gene in self.pathways :
                            continue;
                        else:
                            pway = self.add_pathway_prefix + gene
                            self.pathways[pway] = [gene]
                            if not 'None' in str(type(self.gene_mapping)) and gene in self.gene_mapping.keys():
                                self.pathway_names[pway] = self.gene_mapping[gene]
                            else:
                                self.pathway_names[pway] = gene
                    else :
                        if not self.pathway_tag is None:
                            if not self.pathway_tag in pathway:
                                continue
                        if not 'None' in str(type(self.gene_mapping)) and gene in self.gene_mapping.keys():
                            gene = self.gene_mapping[gene]
                        if pathway in self.pathways:
                            self.pathways[pathway].append(gene)
                        else :
                            if True :
                                self.pathways[pathway] = [gene]
                                self.pathway_names[pathway] = pathway_name
                            else :
                                if len(set(self.lexical_restrictions)-set(lex_restrict.split()))<len(self.lexical_restrictions):
                                    if len( set(self.lexical_restrictions) & set(lex_restrict.split()) )>0:
                                        self.pathways[pathway] = [gene]
                                        self.pathway_names[pathway] = pathway_name

    def add_pathway_synonyms( self , synonym_dict, prefix='' ) :
        for pathway, pathway_name, genes in self.get_generator() :
            synonyms = []
            for g in genes:
                if g in synonym_dict or prefix+g in synonym_dict :
                    k = prefix * ( not g in synonym_dict ) + g
                    sdg = synonym_dict[ k ]
                    if 'list' in str(type(sdg)) :
                        [ synonyms.append(s) for s in sdg ]
                    if 'str' in str(type(sdg)) :
                        synonyms.append(sdg)
            [ self.pathways[pathway].append(s) for s in synonyms ]


    def read_extra_information ( self , file_name , valpos = 0 , keypos = 1 ,
                                 add_value_prefix = '' , required_content = 'HSA' ,
                                 add_desc_pos = None , comp_pos = None ) :

        with open ( file_name ) as input :
            pathway, genes = "", []
            for line in input :
                lspl = line.split( '\t' )
                pw = lspl[ keypos ]
                if required_content in pw :
                    gn = add_value_prefix + lspl[ valpos ]
                    if not self.restrict_id_to is None : 
                        if not pw in self.restrict_id_to :
                            continue
                    if pw in self.pathways :
                        self.pathways[pw] .append(gn)
                    else :
                        self.pathways[pw] = [gn]
                        self.pathway_names[ pw ]  = ''
                    if not add_desc_pos is None :
                        if lspl[add_desc_pos] not in self.pathway_names[pw] :
                            self.pathway_names[pw] += lspl[add_desc_pos] + ', '
                    if not comp_pos is None :
                        if self.pathway_compartments is None :
                            self.pathway_compartments = {}
                        if pw not in self.pathway_compartments :
                            self.pathway_compartments[pw] = [ ''.join(find_compartment(lspl[comp_pos])) ]
                        else :
                            self.pathway_compartments[pw] .append( ''.join(find_compartment(lspl[comp_pos])) )

    def add_extra_information_from_dictionary( self , dictionary, map_onto_genes=True ):
        for pathway in self.pathways :
            genes = self.pathways[pathway]
            add_these = []
            [ [ add_these.append(v) for v in dictionary[g] ] for g in genes if g in dictionary.keys() ]
            if len(add_these)>0:
                [ self.pathways[pathway].append(a) for a in add_these]

def print_generator( generator , show_max_nr=3 , find_str=None ):
    for pathway,pathway_name,genes in generator:
        if not find_str is None:
            if not find_str in pathway:
                continue
        print('Pathway: ' , pathway )
        print('Pathway name: ', pathway_name )
        print('Gene amount : ', len(genes), ' \t Gene name of first: ', genes[0] )
        print('Gene dump   : ', genes )
        show_max_nr-=1
        if show_max_nr == 0:
            break;

def create_listdictionary_from_file( filename , delimiter = '\t' ,add_key_prefix='' ):
    wanted_dictionary = dict()
    with open( filename ) as input :
        value, keyval, descriptor = "", "", []
        for line in input:
            all_vals = line.replace('\n','').replace(' ','').split( delimiter )
            keyv = add_key_prefix+all_vals[0]
            valv = all_vals[1:]
            wanted_dictionary[keyv] = valv
    return ( wanted_dictionary )

def flatten_generator(pathway_generator_object):
    for pathway,genes in pathway_generator_object.pathways.items() :
        ngenes = []
        [ ngenes.append(g) if 'str' in str(type(g)) else [ ngenes.append(q) for q in g ] for g in genes ]
        pathway_generator_object.pathways[pathway] = list(set(ngenes))
        if pathway in pathway_generator_object.pathway_compartments:
            pathway_generator_object.pathway_names[pathway] = '[' + \
                       ','.join(list(set( pathway_generator_object.pathway_compartments[pathway] )) ) + \
                       '] ' + pathway_generator_object.pathway_names[pathway]


def make_compartment_gmt ( reactome_dir ='./', 
        rfiles = ['Ensembl2Reactome_All_Levels.txt','Ensembl2Reactome_PE_All_Levels.txt'] 
    ):
    reac_obj = Reactome ( reactome_dir + rfiles[0] )
    reac_obj .read_extra_information ( reactome_dir + rfiles[1] , add_desc_pos=5, comp_pos=2, keypos=3 )
    flatten_generator ( reac_obj )
    reac_obj .make_gmt_pathway_file('./reactome_pathways_with_compartments.gmt')
    make_group_analytes_unique('./reactome_pathways_with_compartments.gmt')
    return ( reac_obj )

if __name__ == '__main__' :
    print_generator (  make_compartment_gmt(reactome_dir = '/home/richard/Projects/FyFaProjects/ReactomeDB/data/reactome_2019may/').get_generator()  )

    compartment_dict = {}
    if True :
        with open ( './reactome_pathways_with_compartments.gmt' ) as input :
            for line_ in input :
                if not '[' in line_ or not ']' in line_:
                    continue
                line = line_.replace('\n','') ; lspl = line.split('\t')
                analytes = lspl[2:]
                all_comp_entries = lspl[1].split('[')[-1].split(']')[0].split(',')
                for comp in set(all_comp_entries) :
                    if comp in compartment_dict :
                        oalas = compartment_dict[comp]
                        [ oalas.append(c) for c in set(analytes) ]
                        compartment_dict[ comp ] = list(set(analytes))
                    else :
                        compartment_dict[ comp ] = list(set(analytes))
    if True :
        f_i = './Ensembl2Reactome_PE_All_Levels.txt'
        fname_o = 'compartment_genes.gmt'

        with open( f_i ) as input:
            for line in input:
                if 'ENSG' in line and 'HSA' in line:
                    lsp = line.split('\t')
                    g = lsp[0]; comp = lsp[2].split(']')[0].split('[')[1]
                    if not comp in compartment_dict :
                        compartment_dict[comp] = [g]
                    else :
                        g_ = compartment_dict[comp] ; g_ .append( g )
                        compartment_dict[comp] = list( set(g_) )
        for k,v in compartment_dict.items() :
            if True :
                v_ = [ w for w in v if 'ENSG' in w ]
            print ( k+'\t'+k+' from Ensembl2Reactome_PE_All_Levels\t'+'\t'.join(v_) , file = open ( fname_o,'a' ) )

