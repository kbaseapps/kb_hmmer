#BEGIN_HEADER
import os
import sys
import shutil
import hashlib
import subprocess
import requests
import re
import traceback
import uuid
from datetime import datetime
from pprint import pprint, pformat
import numpy as np
import math
import gzip

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from biokbase.workspace.client import Workspace as workspaceService
from requests_toolbelt import MultipartEncoder
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

#END_HEADER


class kb_hmmer:
    '''
    Module Name:
    kb_hmmer

    Module Description:
    ** A KBase module: kb_hmmer
**
** This module contains HMMER Hidden Markov Model Sequence Search and Alignment
**
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    workspaceURL = None
    shockURL = None
    handleURL = None

    HMMER_BUILD      = '/kb/module/hmmer/binaries/hmmbuild'  # construct profile HMM(s) from MSA(s)
    HMMER_MAKE_DB    = '/kb/module/hmmer/binaries/makehmmerdb'  # build a HMMER binary db from a seq file
    HMMER_SEARCH     = '/kb/module/hmmer/binaries/hmmsearch'  # search profile(s) against a sequence db

    HMMER_PHMMER     = '/kb/module/hmmer/binaries/phmmer'  # search protein sequence(s) against a protein sequence db
    HMMER_NHMMER     = '/kb/module/hmmer/binaries/nhmmer'  # search nuc sequence(s) against a nuc sequence db
    HMMER_JACKHAMMER = '/kb/module/hmmer/binaries/jackhmmer'  # iteratively search sequence(s) against a protein db

    #HMMER_ALIGN      = '/kb/module/hmmer/binaries/hmmalign'  # align sequences to a profile HMM
    #HMMER_PRESS      = '/kb/module/hmmer/binaries/hmmpress'  # prepare HMM db for hmmscan
    #HMMER_SCAN       = '/kb/module/hmmer/binaries/hmmscan'  # scan prot sequence(s) against protein profile db
    #HMMER_NSCAN       = '/kb/module/hmmer/binaries/nhmmscan'  # scan nuc sequence(s) against nuc profile db


    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    def get_single_end_read_library(self, ws_data, ws_info, forward):
        pass

    def get_feature_set_seqs(self, ws_data, ws_info):
        pass

    def get_genome_feature_seqs(self, ws_data, ws_info):
        pass

    def get_genome_set_feature_seqs(self, ws_data, ws_info):
        pass

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.scratch = os.path.abspath(config['scratch'])
        # HACK!! temporary hack for issue where megahit fails on mac because of silent named pipe error
        #self.host_scratch = self.scratch
        self.scratch = os.path.join('/kb','module','local_scratch')
        # end hack
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        #END_CONSTRUCTOR
        pass


    # Helper script borrowed from the transform service, logger removed
    #
    def upload_file_to_shock(self,
                             console,  # DEBUG
                             shock_service_url = None,
                             filePath = None,
                             ssl_verify = True,
                             token = None):
        """
        Use HTTP multi-part POST to save a file to a SHOCK instance.
        """
        self.log(console,"UPLOADING FILE "+filePath+" TO SHOCK")

        if token is None:
            raise Exception("Authentication token required!")

        #build the header
        header = dict()
        header["Authorization"] = "Oauth {0}".format(token)
        if filePath is None:
            raise Exception("No file given for upload to SHOCK!")

        dataFile = open(os.path.abspath(filePath), 'rb')
        m = MultipartEncoder(fields={'upload': (os.path.split(filePath)[-1], dataFile)})
        header['Content-Type'] = m.content_type

        #logger.info("Sending {0} to {1}".format(filePath,shock_service_url))
        try:
            response = requests.post(shock_service_url + "/node", headers=header, data=m, allow_redirects=True, verify=ssl_verify)
            dataFile.close()
        except:
            dataFile.close()
            raise
        if not response.ok:
            response.raise_for_status()
        result = response.json()
        if result['error']:
            raise Exception(result['error'][0])
        else:
            return result["data"]


    def upload_SingleEndLibrary_to_shock_and_ws (self,
                                                 ctx,
                                                 console,  # DEBUG
                                                 workspace_name,
                                                 obj_name,
                                                 file_path,
                                                 provenance,
                                                 sequencing_tech):

        self.log(console,'UPLOADING FILE '+file_path+' TO '+workspace_name+'/'+obj_name)

        # 1) upload files to shock
        token = ctx['token']
        forward_shock_file = self.upload_file_to_shock(
            console,  # DEBUG
            shock_service_url = self.shockURL,
            filePath = file_path,
            token = token
            )
        #pprint(forward_shock_file)
        self.log(console,'SHOCK UPLOAD DONE')

        # 2) create handle
        self.log(console,'GETTING HANDLE')
        hs = HandleService(url=self.handleURL, token=token)
        forward_handle = hs.persist_handle({
                                        'id' : forward_shock_file['id'], 
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': forward_shock_file['file']['name'],
                                        'remote_md5': forward_shock_file['file']['checksum']['md5']})

        
        # 3) save to WS
        self.log(console,'SAVING TO WORKSPACE')
        single_end_library = {
            'lib': {
                'file': {
                    'hid':forward_handle,
                    'file_name': forward_shock_file['file']['name'],
                    'id': forward_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':forward_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fasta',
                'size':forward_shock_file['file']['size']
            },
            'sequencing_tech':sequencing_tech
        }
        self.log(console,'GETTING WORKSPACE SERVICE OBJECT')
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        self.log(console,'SAVE OPERATION...')
        new_obj_info = ws.save_objects({
                        'workspace':workspace_name,
                        'objects':[
                            {
                                'type':'KBaseFile.SingleEndLibrary',
                                'data':single_end_library,
                                'name':obj_name,
                                'meta':{},
                                'provenance':provenance
                            }]
                        })
        self.log(console,'SAVED TO WORKSPACE')

        return new_obj_info[0]

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.scratch = os.path.abspath(config['scratch'])
        # HACK!! temporary hack for issue where megahit fails on mac because of silent named pipe error
        #self.host_scratch = self.scratch
        self.scratch = os.path.join('/kb','module','local_scratch')
        # end hack
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        #END_CONSTRUCTOR
        pass


    def HMMER_MSA_Search(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN HMMER_MSA_Search
        console = []
        self.log(console,'Running HMMER_MSA_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running HMMER_MSA_Search with params='
#        report += "\n"+pformat(params)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'input_msa_name' not in params:
            raise ValueError('input_msa_name parameter is required')
        if 'input_many_name' not in params:
            raise ValueError('input_many_name parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')


        #### Get the input_msa object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_msa_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            input_msa_type = info[2].split('.')[1].split('-')[0]

        except Exception as e:
            raise ValueError('Unable to fetch input_msa_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        if input_msa_type == 'MSA':
            MSA_in = data
            row_order = []
            default_row_labels = dict()
            if 'row_order' in MSA_in.keys():
                row_order = MSA_in['row_order']
            else:
                row_order = sorted(MSA_in['alignment'].keys())

            if 'default_row_labels' in MSA_in.keys():
                default_row_labels = MSA_in['default_row_labels']
            else:
                for row_id in row_order:
                    default_row_labels[row_id] = row_id

            # export features to CLUSTAL formatted MSA (HMMER BUILD seems to only take CLUSTAL)
            input_MSA_file_path = os.path.join(self.scratch, params['input_msa_name']+".fasta")
            self.log(console, 'writing MSA file: '+input_MSA_file_path)

            # set header
            header = 'CLUSTAL W (1.81) multiple sequence alignment'

            # get longest id
            longest_row_id_len = 0
            for row_id in row_order:
                if len(row_id) > longest_row_id_len:
                    longest_row_id_len = len(row_id)
            # make sure rows are all same length
            row_id_0 = row_order[0]
            row_len = len(MSA_in['alignment'][row_id_0])
            for row_id in row_order:
                if len(MSA_in['alignment'][row_id]) != row_len:
                    raise ValueError("MSA alignment rows are not constant length")
            # get alignment line (just storing identity markers)
            conservation_symbol = []
            for i in range(row_len):
                first_seen_char = MSA_in['alignment'][row_id_0][i]
                symbol = '*'
                for row_id in row_order:
                    if MSA_in['alignment'][row_id][i] == '-' or MSA_in['alignment'][row_id][i] != first_seen_char:
                        symbol = ' '
                        break
                conservation_symbol.append(symbol)

            # break up MSA into 60 char chunks
            records = []
            whole_chunks = int(math.floor(row_len/60))
            if whole_chunks > 0:
                for j in range(whole_chunks):
                    for row_id in row_order:
                        padding = ''
                        if longest_row_id_len-len(row_id) > 0:
                            for i in range(0,longest_row_id_len-len(row_id)):
                                padding += ' '
                        records.append(row_id + padding + " " +
                                       MSA_in['alignment'][row_id][j*60:(j+1)*60])
                    records.append(''.join([' ' for s in range(longest_row_id_len)]) + " " +
                                   conservation_symbol[j*60:(j+1)*60])
            # add incomplete rows
            j=whole_chunks
            for row_id in row_order:
                padding = ''
                if longest_row_id_len-len(row_id) > 0:
                    for i in range(0,longest_row_id_len-len(row_id)):
                        padding += ' '
                records.append(row_id + padding + " " +
                               MSA_in['alignment'][row_id][j*60:row_len])
                records.append(''.join([' ' for s in range(longest_row_id_len)]) + " " +
                               conservation_symbol[j*60:row_len])
            
            # write that sucker
            with open(input_MSA_file_path,'w',0) as input_MSA_file_handle:
                input_MSA_file_handle.write(header+"\n\n")
                input_MSA_file_handle.write("\n".join(records)+"\n")

            # DEBUG
            report += "MSA:\n"
            report += header+"\n\n"
            report += "\n".join(records)+"\n"
            self.log(console,report)


            # Determine whether nuc or protein sequences
            #
            NUC_MSA_pattern = re.compile("^[\.\-_ACGTUXNRYSWKMBDHVacgtuxnryswkmbdhv \t\n]+$")
            all_seqs_nuc = True
            for row_id in row_order:
                #self.log(console, row_id+": '"+MSA_in['alignment'][row_id]+"'")
                if NUC_MSA_pattern.match(MSA_in['alignment'][row_id]) == None:
                    all_seqs_nuc = False
                    break

        # Missing proper input_type
        #
        else:
            raise ValueError('Cannot yet handle input_name type of: '+type_name)

        #### Get the input_many object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_many_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            many_type_name = info[2].split('.')[1].split('-')[0]

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        # Handle overloading (input_many can be FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_many_featureSet = data

            genome2Features = {}
            features = input_many_featureSet['elements']
            for fId in features.keys():
                genomeRef = features[fId][0]
                if genomeRef not in genome2Features:
                    genome2Features[genomeRef] = []
                genome2Features[genomeRef].append(fId)

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)
            records = []
            feature_written = dict()
            for genomeRef in genome2Features:
                genome = ws.get_objects([{'ref':genomeRef}])[0]['data']
                these_genomeFeatureIds = genome2Features[genomeRef]
                for feature in genome['features']:
                    if feature['id'] in these_genomeFeatureIds:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG

                            # HMMER SEARCH is prot-prot in this implementation
                            if feature['type'] != 'CDS':
                                self.log(console,"skipping non-CDS feature "+feature['id'])
                                continue
                            elif 'protein_translation' not in feature or feature['protein_translation'] == None:
                                self.log(console,"bad CDS feature "+feature['id'])
                                raise ValueError("bad CDS feature "+feature['id'])
                            else:
                                #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                                record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                                records.append(record)
            SeqIO.write(records, many_forward_reads_file_path, "fasta")


        # Genome
        #
        elif many_type_name == 'Genome':
            input_many_genome = data
            input_many_genome_ref = str(info[6])+'/'+str(info[0])+'/'+str(info[4])

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)
            records = []
            feature_written = dict()
            for feature in input_many_genome['features']:
                try:
                    f_written = feature_written[feature['id']]
                except:
                    feature_written[feature['id']] = True
                    #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG

                    # HMMER SEARCH is prot-prot in this implementation
                    #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=input_many_genome['id'])
                    if feature['type'] != 'CDS':
                        #self.log(console,"skipping non-CDS feature "+feature['id'])  # too much chatter for a Genome
                        continue
                    elif 'protein_translation' not in feature or feature['protein_translation'] == None:
                        self.log(console,"bad CDS feature "+feature['id'])
                        raise ValueError("bad CDS feature "+feature['id'])
                    else:
                        record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=input_many_genome['id'])
                        records.append(record)
            SeqIO.write(records, many_forward_reads_file_path, "fasta")


        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = data

            # export features to FASTA file
            many_forward_reads_file_path = os.path.join(self.scratch, params['input_many_name']+".fasta")
            self.log(console, 'writing fasta file: '+many_forward_reads_file_path)

            records = []
            feature_written = dict()
            for genome_name in input_many_genomeSet['elements'].keys():
                if 'ref' in input_many_genomeSet['elements'][genome_name] and \
                         input_many_genomeSet['elements'][genome_name]['ref'] != None:
                    genome = ws.get_objects([{'ref': input_many_genomeSet['elements'][genome_name]['ref']}])[0]['data']
                    for feature in genome['features']:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # HMMER SEARCH is prot-prot in this implementation
                            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            if feature['type'] != 'CDS':
                                #self.log(console,"skipping non-CDS feature "+feature['id'])  # too much chatter for a Genome
                                continue
                            elif 'protein_translation' not in feature or feature['protein_translation'] == None:
                                self.log(console,"bad CDS feature "+feature['id'])
                                raise ValueError("bad CDS feature "+feature['id'])
                            else:
                                record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                                records.append(record)

                elif 'data' in input_many_genomeSet['elements'][genome_name] and \
                        input_many_genomeSet['elements'][genome_name]['data'] != None:
                    genome = input_many_genomeSet['elements'][genome_name]['data']
                    for feature in genome['features']:
                        try:
                            f_written = feature_written[feature['id']]
                        except:
                            feature_written[feature['id']] = True
                            #self.log(console,"kbase_id: '"+feature['id']+"'")  # DEBUG
                            # HMMER SEARCH is prot-prot in this implementation
                            #record = SeqRecord(Seq(feature['dna_sequence']), id=feature['id'], description=genome['id'])
                            if feature['type'] != 'CDS':
                                continue
                            elif 'protein_translation' not in feature or feature['protein_translation'] == None:
                                self.log(console,"bad CDS feature "+feature['id'])
                                raise ValueError("bad CDS feature "+feature['id'])
                            else:
                                record = SeqRecord(Seq(feature['protein_translation']), id=feature['id'], description=genome['id'])
                                records.append(record)

                else:
                    raise ValueError('genome '+genome_name+' missing')

            SeqIO.write(records, many_forward_reads_file_path, "fasta")
            
        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)            


        # Build HMM from MSA
        #
        # SYNTAX (from http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf)
        #
        # hmmbuild --informat fasta <hmmfile.out> <msafile>
        #
        hmmer_build_bin = self.HMMER_BUILD
        hmmer_build_cmd = [hmmer_build_bin]

        # check for necessary files
        if not os.path.isfile(hmmer_build_bin):
            raise ValueError("no such file '"+hmmer_build_bin+"'")
        if not os.path.isfile(input_MSA_file_path):
            raise ValueError("no such file '"+input_MSA_file_path+"'")
        elif not os.path.getsize(input_MSA_file_path) > 0:
            raise ValueError("empty file '"+input_MSA_file_path+"'")

        HMM_file_path = input_MSA_file_path+".HMM"

        hmmer_build_cmd.append('--informat')
        hmmer_build_cmd.append('FASTA')
        hmmer_build_cmd.append(HMM_file_path)
        hmmer_build_cmd.append(input_MSA_file_path)

        # Run HMMER_BUILD, capture output as it happens
        #
        self.log(console, 'RUNNING HMMER_BUILD:')
        self.log(console, '    '+' '.join(hmmer_build_cmd))
#        report += "\n"+'running HMMER_BUILD:'+"\n"
#        report += '    '+' '.join(hmmer_build_cmd)+"\n"

        p = subprocess.Popen(hmmer_build_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

        while True:
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running HMMER_BUILD, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))

        # Check for HMM output
        if not os.path.isfile(HMM_file_path):
            raise ValueError("HMMER_BUILD failed to create HMM file '"+HMM_file_path+"'")
        elif not os.path.getsize(HMM_file_path+"'") > 0:
            raise ValueError("HMMER_BUILD created empty HMM file '"+HMM_file_path+"'")


        ### Construct the HMMER_SEARCH command
        #
        # SYNTAX (from http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf)
        #
        # hmmsearch --tblout -o <output> -A <MSA_out> --noali --notextw -E <e_value> -T <bit_score> <hmmfile> <seqdb>
        #
        hmmer_search_bin = self.HMMER_SEARCH
        hmmer_search_cmd = [hmmer_search_bin]

        # check for necessary files
        if not os.path.isfile(hmmer_search_bin):
            raise ValueError("no such file '"+hmmer_search_bin+"'")
        if not os.path.isfile(HMM_file_path):
            raise ValueError("no such file '"+HMM_file_path+"'")
        elif not os.path.getsize(HMM_file_path):
            raise ValueError("empty file '"+HMM_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")
        elif not os.path.getsize(many_forward_reads_file_path):
            raise ValueError("empty file '"+many_forward_reads_file_path+"'")

        # set the output path
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_hit_TAB_file_path = os.path.join(output_dir, 'hitout.txt');
        output_hit_MSA_file_path = os.path.join(output_dir, 'msaout.txt');
        output_filtered_fasta_file_path = os.path.join(output_dir, 'output_filtered.fasta');

        # this is command for basic search mode
        hmmer_search_cmd.append('--tblout')
        hmmer_search_cmd.append('-o')
        hmmer_search_cmd.append(output_TAB_file_path)
        hmmer_search_cmd.append('-A')
        hmmer_search_cmd.append(output_MSA_file_path)
        hmmer_search_cmd.append('--noali')
        hmmer_search_cmd.append('--notextw')
        hmmer_search_cmd.append('-E')
        hmmer_search_cmd.append(str(params['e_value']))
        hmmer_search_cmd.append(HMM_file_path)
        hmmer_search_cmd.append(many_forward_reads_file_path)

        # options
#        if 'maxaccepts' in params:
#            if params['maxaccepts']:
#                hmmer_search_cmd.append('-max_target_seqs')
#                hmmer_search_cmd.append(str(params['maxaccepts']))

        # Run HMMER, capture output as it happens
        #
        self.log(console, 'RUNNING HMMER_SEARCH:')
        self.log(console, '    '+' '.join(hmmer_search_cmd))
#        report += "\n"+'running HMMER_SEARCH:'+"\n"
#        report += '    '+' '.join(hmmer_search_cmd)+"\n"

        p = subprocess.Popen(hmmer_search_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

        while True:
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running HMMER_SEARCH, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))


        # Check for output
        if not os.path.isfile(output_hit_TAB_file_path):
            raise ValueError("HMMER_SEARCH failed to create TAB file '"+output_hit_TAB_file_path+"'")
        elif not os.path.getsize(output_hit_TAB_file_path+"'") > 0:
            raise ValueError("HMMER_SEARCH created empty TAB file '"+output_hit_TAB_file_path+"'")
        if not os.path.isfile(output_hit_MSA_file_path):
            raise ValueError("HMMER_SEARCH failed to create MSA file '"+output_hit_MSA_file_path+"'")
        elif not os.path.getsize(output_hit_MSA_file_path+"'") > 0:
            raise ValueError("HMMER_SEARCH created empty MSA file '"+output_hit_MSA_file_path+"'")


        # DEBUG
        report = "TAB:\n\n"
        with open (output_hit_TAB_file_path, 'r') as output_handle:
            for line in output_handle:
                report += line+"\n"
        report += "\n\nMSA:\n\n"
        with open (output_hit_MSA_file_path, 'r') as output_handle:
            for line in output_handle:
                report += line+"\n"


        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG

        reportObj = {
            'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':'HMMER_MSA_Search hits'}],
            'text_message':report
        }

        reportName = 'hmmer_report_'+str(hex(uuid.getnode()))
        report_obj_info = ws.save_objects({
#                'id':info[6],
                'workspace':params['workspace_name'],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        #'provenance':provenance  # DEBUG
                        'provenance':None
                    }
                ]
            })[0]

        self.log(console,"BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
        self.log(console,"HMMER_MSA_Search DONE")
        #END HMMER_MSA_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method HMMER_MSA_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
