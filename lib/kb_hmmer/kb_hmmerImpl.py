# -*- coding: utf-8 -*-
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
#from biokbase.workspace.client import Workspace as workspaceService
from Workspace.WorkspaceClient import Workspace as workspaceService
from requests_toolbelt import MultipartEncoder
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

# SDK Utils
from KBaseDataObjectToFileUtils.KBaseDataObjectToFileUtilsClient import KBaseDataObjectToFileUtils
from DataFileUtil.DataFileUtilClient import DataFileUtil as DFUClient
from KBaseReport.KBaseReportClient import KBaseReport

# silence whining
import requests
requests.packages.urllib3.disable_warnings()

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

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.1.2"
    GIT_URL = "https://github.com/kbaseapps/kb_hmmer"
    GIT_COMMIT_HASH = "8c8caa67f5a699fefb9a1f4f8590232d0b779ef0"

    #BEGIN_CLASS_HEADER
    workspaceURL = None
    shockURL     = None
    handleURL    = None
    callbackURL  = None
    scratch      = None

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
                        })[0]
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
        self.serviceWizardURL = config['service-wizard-url']

#        self.callbackURL = os.environ['SDK_CALLBACK_URL'] if os.environ['SDK_CALLBACK_URL'] != None else 'https://kbase.us/services/njs_wrapper'
        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
        if self.callbackURL == None:
            raise ValueError ("SDK_CALLBACK_URL not set in environment")

        self.scratch = os.path.abspath(config['scratch'])
        if self.scratch == None:
            self.scratch = os.path.join('/kb','module','local_scratch')
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        # set i/o dirs
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        self.input_dir = os.path.join(self.scratch,'input.'+str(timestamp))
        self.output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(self.input_dir):
            os.makedirs(self.input_dir)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        #END_CONSTRUCTOR
        pass


    def HMMER_MSA_Search(self, ctx, params):
        """
        Methods for HMMER search of an MSA against many sequences 
        **
        **    overloading as follows:
        **        input_msa_ref: MSA
        **        input_many_ref: SingleEndLibrary, FeatureSet, Genome, GenomeSet
        **        output_name: SingleEndLibrary (if input_many is SELib), (else) FeatureSet
        :param params: instance of type "HMMER_Params" (HMMER Input Params)
           -> structure: parameter "workspace_name" of type "workspace_name"
           (** The workspace object refs are of form: ** **    objects =
           ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_many_ref" of type "data_obj_ref", parameter
           "input_msa_ref" of type "data_obj_ref", parameter
           "output_filtered_name" of type "data_obj_name", parameter
           "e_value" of Double, parameter "bitscore" of Double, parameter
           "maxaccepts" of Double
        :returns: instance of type "HMMER_Output" (HMMER Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN HMMER_MSA_Search
        console = []
        invalid_msgs = []
        search_tool_name = 'HMMER_MSA_prot'
        self.log(console,'Running '+search_tool_name+'_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running '+search_tool_name+'_Search with params='
#        report += "\n"+pformat(params)
        #appropriate_sequence_found_in_one_input = False
        appropriate_sequence_found_in_MSA_input = False
        appropriate_sequence_found_in_many_input = False
        genome_id_feature_id_delim = '.f:'

        # set hmmer_dir
        hmmer_dir = os.path.join(self.output_dir, 'hmmer_run')
        if not os.path.exists(hmmer_dir):
            os.makedirs(hmmer_dir)


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
#        if 'input_one_ref' not in params:
#            raise ValueError('input_one_ref parameter is required')
        if 'input_msa_ref' not in params:
            raise ValueError('input_msa_ref parameter is required')
        if 'input_many_ref' not in params:
            raise ValueError('input_many_ref parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')


        # set local names
#        input_one_ref = params['input_one_ref']
        input_msa_ref = params['input_msa_ref']
        input_many_ref = params['input_many_ref']
        

        #### Get the input_msa object
        ##
#        if input_one_feature_id == None:
#            self.log(invalid_msgs,"input_one_feature_id was not obtained from Query Object: "+input_one_name)
#        master_row_idx = 0
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            #objects = ws.get_objects([{'ref': input_msa_ref}])
            objects = ws.get_objects2({'objects':[{'ref': input_msa_ref}]})['data']
            input_msa_data = objects[0]['data']
            info = objects[0]['info']
            input_msa_name = str(info[1])
            msa_type_name = info[2].split('.')[1].split('-')[0]

        except Exception as e:
            raise ValueError('Unable to fetch input_msa_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        if msa_type_name == 'MSA':
            MSA_in = input_msa_data
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

            # determine row index of query sequence
#            for row_id in row_order:
#                master_row_idx += 1
#                if row_id == input_one_feature_id:
#                    break
#            if master_row_idx == 0:
#                self.log(invalid_msgs,"Failed to find query id "+input_one_feature_id+" from Query Object "+input_one_name+" within MSA: "+input_msa_name)

            
            # export features to CLUSTAL formatted MSA (HMMER BUILD seems to only take CLUSTAL)
            input_MSA_file_path = os.path.join(hmmer_dir, input_msa_name+".clustal")
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
            conservation_symbol = ''
            for i in range(row_len):
                first_seen_char = MSA_in['alignment'][row_id_0][i]
                symbol = '*'
                for row_id in row_order:
                    if MSA_in['alignment'][row_id][i] == '-' or MSA_in['alignment'][row_id][i] != first_seen_char:
                        symbol = ' '
                        break
                conservation_symbol += symbol

            # break up MSA into 60 char chunks
            records = []
            chunk_len = 60
            whole_chunks = int(math.floor(row_len/chunk_len))
            if whole_chunks > 0:
                for j in range(whole_chunks):
                    records.append('')
                    for row_id in row_order:
                        padding = ''
                        if longest_row_id_len-len(row_id) > 0:
                            for i in range(0,longest_row_id_len-len(row_id)):
                                padding += ' '
                        records.append(row_id + padding + " " +
                                       MSA_in['alignment'][row_id][j*chunk_len:(j+1)*chunk_len])
                    records.append(''.join([' ' for s in range(longest_row_id_len)]) + " " +
                                   conservation_symbol[j*chunk_len:(j+1)*chunk_len])

            # add final rows
            if (row_len % chunk_len) != 0:
                j=whole_chunks
                records.append('')
                for row_id in row_order:
                    padding = ''
                    if longest_row_id_len-len(row_id) > 0:
                        for i in range(0,longest_row_id_len-len(row_id)):
                            padding += ' '
                    records.append(row_id + padding + " " +
                                   MSA_in['alignment'][row_id][j*chunk_len:row_len])
                records.append(''.join([' ' for s in range(longest_row_id_len)]) + " " +
                               conservation_symbol[j*chunk_len:row_len])
            
            # write that sucker
            with open(input_MSA_file_path,'w',0) as input_MSA_file_handle:
                input_MSA_file_handle.write(header+"\n")
                input_MSA_file_handle.write("\n".join(records)+"\n")

            # DEBUG
            #report += "MSA:\n"
            #report += header+"\n"
            #report += "\n".join(records)+"\n"
            #self.log(console,report)


            # Determine whether nuc or protein sequences
            #
            self.log (console, "CHECKING MSA for PROTEIN seqs...")  # DEBUG
            PROT_MSA_pattern = re.compile("^[\.\-_acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
            #NUC_MSA_pattern = re.compile("^[\.\-_ACGTUXNRYSWKMBDHVacgtuxnryswkmbdhv \t\n]+$")
            appropriate_sequence_found_in_MSA_input = True
            for row_id in row_order:
                #self.log(console, row_id+": '"+MSA_in['alignment'][row_id]+"'")    # DEBUG
                if not PROT_MSA_pattern.match(MSA_in['alignment'][row_id]):
                    self.log(invalid_msgs,"BAD record for MSA row_id: "+row_id+"\n"+MSA_in['alignment'][row_id]+"\n")
                    appropriate_sequence_found_in_MSA_input = False
                    break

        # Missing proper input_type
        #
        else:
            raise ValueError('Cannot yet handle input_msa type of: '+msa_type_name)


        #### Get the input_many object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            #objects = ws.get_objects([{'ref': input_many_ref}])
            objects = ws.get_objects2({'objects':[{'ref': input_many_ref}]})['data']
            input_many_data = objects[0]['data']
            info = objects[0]['info']
            input_many_name = str(info[1])
            many_type_name = info[2].split('.')[1].split('-')[0]

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()


        # Handle overloading (input_many can be SequenceSet, FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'SequenceSet':
            try:
                input_many_sequenceSet = input_many_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get SequenceSet: ' + str(e))

            header_id = input_many_sequenceSet['sequences'][0]['sequence_id']
            many_forward_reads_file_path = os.path.join(self.output_dir, header_id+'.fasta')
            many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(many_forward_reads_file_path))

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                sequence_str = seq_obj['sequence']

                PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
                #DNA_pattern = re.compile("^[acgtuACGTUnryNRY ]+$")
                if not PROT_pattern.match(sequence_str):
                    self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
                    continue
                appropriate_sequence_found_in_many_input = True
                many_forward_reads_file_handle.write('>'+header_id+"\n")
                many_forward_reads_file_handle.write(sequence_str+"\n")
            many_forward_reads_file_handle.close();
            self.log(console, 'done')


        # FeatureSet
        #
        if many_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_many_featureSet = input_many_data
            many_forward_reads_file_dir = self.output_dir
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            FeatureSetToFASTA_params = {
                'featureSet_ref':      input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            #SERVICE_VER = 'release'
            SERVICE_VER = 'dev'
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            many_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            feature_ids_by_genome_ref = FeatureSetToFASTA_retVal['feature_ids_by_genome_ref']
            if len(feature_ids_by_genome_ref.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Genome
        #
        elif many_type_name == 'Genome':
            many_forward_reads_file_dir = self.output_dir
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeToFASTA_params = {
                'genome_ref':          input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            #SERVICE_VER = 'release'
            SERVICE_VER = 'dev'
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
            GenomeToFASTA_retVal = DOTFU.GenomeToFASTA (GenomeToFASTA_params)
            many_forward_reads_file_path = GenomeToFASTA_retVal['fasta_file_path']
            feature_ids = GenomeToFASTA_retVal['feature_ids']
            if len(feature_ids) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "Genome2Fasta() took "+str(end_time-beg_time)+" secs")
            

        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = input_many_data
            many_forward_reads_file_dir = self.output_dir
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeSetToFASTA_params = {
                'genomeSet_ref':       input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            #SERVICE_VER = 'release'
            SERVICE_VER = 'dev'
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
            GenomeSetToFASTA_retVal = DOTFU.GenomeSetToFASTA (GenomeSetToFASTA_params)
            many_forward_reads_file_path = GenomeSetToFASTA_retVal['fasta_file_path_list'][0]
            feature_ids_by_genome_id = GenomeSetToFASTA_retVal['feature_ids_by_genome_id']
            if len(feature_ids_by_genome_id.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+many_type_name)            

        # check for failed input file creation
        #
#        if not appropriate_sequence_found_in_one_input:
#            self.log(invalid_msgs,"no protein sequences found in '"+input_one_name+"'")
        if not appropriate_sequence_found_in_MSA_input:
            self.log(invalid_msgs,"no protein sequences found in '"+input_msa_name+"'")
        if not appropriate_sequence_found_in_many_input:
            self.log(invalid_msgs,"no protein sequences found in '"+input_many_name+"'")


        # input data failed validation.  Need to return
        #
        if len(invalid_msgs) > 0:

            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
#            provenance[0]['input_ws_objects'].append(input_one_ref)
            provenance[0]['input_ws_objects'].append(input_msa_ref)
            provenance[0]['input_ws_objects'].append(input_many_ref)
            provenance[0]['service'] = 'kb_hmmer'
            provenance[0]['method'] = search_tool_name+'_Search'


            # build output report object
            #
            self.log(console,"BUILDING REPORT")  # DEBUG
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            reportName = 'hmmer_report_'+str(uuid.uuid4())
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            report_obj_info = ws.save_objects({
                    #'id':info[6],
                    'workspace':params['workspace_name'],
                    'objects':[
                        {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance  # DEBUG
                        }
                        ]
                    })[0]

            self.log(console,"BUILDING RETURN OBJECT")
            returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
            self.log(console,search_tool_name+"_Search DONE")
            return [returnVal]


        # Set output paths
        #output_aln_file_path = os.path.join(hmmer_dir, 'alnout.txt');
        #output_extra_file_path = os.path.join(hmmer_dir, 'alnout_extra.txt');
        #output_filtered_fasta_file_path = os.path.join(hmmer_dir, 'output_filtered.faa');


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
        hmmer_build_cmd.append('CLUSTAL')
        hmmer_build_cmd.append(HMM_file_path)
        hmmer_build_cmd.append(input_MSA_file_path)

        # Run HMMER_BUILD, capture output as it happens
        #
        self.log(console, 'RUNNING HMMER_BUILD:')
        self.log(console, '    '+' '.join(hmmer_build_cmd))
#        report += "\n"+'running HMMER_BUILD:'+"\n"
#        report += '    '+' '.join(hmmer_build_cmd)+"\n"

        p = subprocess.Popen(hmmer_build_cmd, \
                             cwd = self.output_dir, \
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
        elif not os.path.getsize(HMM_file_path) > 0:
            raise ValueError("HMMER_BUILD created empty HMM file '"+HMM_file_path+"'")

        # DEBUG
        with open (HMM_file_path, 'r') as HMM_file_handle:
            for line in HMM_file_handle.readlines():
                self.log(console, "HMM_FILE: '"+str(line)+"'")


        ### Construct the HMMER_SEARCH command
        #
        # SYNTAX (from http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf)
        #
        # hmmsearch --tblout <TAB_out> -A <MSA_out> --noali --notextw -E <e_value> -T <bit_score> <hmmfile> <seqdb>
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

        output_hit_TAB_file_path = os.path.join(hmmer_dir, 'hitout.txt');
        output_hit_MSA_file_path = os.path.join(hmmer_dir, 'msaout.txt');
        output_filtered_fasta_file_path = os.path.join(hmmer_dir, 'output_filtered.fasta');

        # this is command for basic search mode
        hmmer_search_cmd.append('--tblout')
        hmmer_search_cmd.append(output_hit_TAB_file_path)
        hmmer_search_cmd.append('-A')
        hmmer_search_cmd.append(output_hit_MSA_file_path)
        hmmer_search_cmd.append('--noali')
        hmmer_search_cmd.append('--notextw')
        hmmer_search_cmd.append('-E')  # can't use -T with -E, so we'll use -E
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
                             cwd = self.output_dir, \
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
        elif not os.path.getsize(output_hit_TAB_file_path) > 0:
            raise ValueError("HMMER_SEARCH created empty TAB file '"+output_hit_TAB_file_path+"'")
        if not os.path.isfile(output_hit_MSA_file_path):
            raise ValueError("HMMER_SEARCH failed to create MSA file '"+output_hit_MSA_file_path+"'")
        elif not os.path.getsize(output_hit_MSA_file_path) > 0:
            raise ValueError("HMMER_SEARCH created empty MSA file '"+output_hit_MSA_file_path+"'")


        # DEBUG
#        report = "TAB:\n\n"
#        with open (output_hit_TAB_file_path, 'r') as output_handle:
#            for line in output_handle:
#                report += line+"\n"
#        report += "\n\nMSA:\n\n"
#        with open (output_hit_MSA_file_path, 'r') as output_handle:
#            for line in output_handle:
#                report += line+"\n"


        # Parse the HMMER tabular output and store ids to filter many set to make filtered object to save back to KBase
        #
        self.log(console, 'PARSING HMMER SEARCH OUTPUT')
        if not os.path.isfile(output_hit_TAB_file_path):
            raise ValueError("failed to create HMMER output: "+output_hit_TAB_file_path)
        elif not os.path.getsize(output_hit_TAB_file_path) > 0:
            raise ValueError("created empty file for HMMER output: "+output_hit_TAB_file_path)
        hit_seq_ids = dict()
        accept_fids = dict()
        output_hit_TAB_file_handle = open (output_hit_TAB_file_path, "r", 0)
        output_aln_buf = output_hit_TAB_file_handle.readlines()
        output_hit_TAB_file_handle.close()
        hit_total = 0
        high_bitscore_line = dict()
        high_bitscore_score = dict()
        high_bitscore_ident = dict()
        high_bitscore_alnlen = dict()
        hit_order = []
        hit_buf = []
        #header_done = False
        for line in output_aln_buf:
            if line.startswith('#'):
                #if not header_done:
                #    hit_buf.append(line)
                continue
            #header_done = True
            #self.log(console,'HIT LINE: '+line)  # DEBUG
            hit_info = re.split ('\s+', line)
            hit_seq_id            = hit_info[0]
            hit_accession         = hit_info[1]
            query_name            = hit_info[2]
            query_accession       = hit_info[3]
            hit_e_value           = float(hit_info[4])
            hit_bitscore          = float(hit_info[5])
            hit_bias              = float(hit_info[6])
            hit_e_value_best_dom  = float(hit_info[7])
            hit_bitscore_best_dom = float(hit_info[8])
            hit_bias_best_dom     = float(hit_info[9])
            hit_expected_dom_n    = float(hit_info[10])
            hit_regions           = float(hit_info[11])
            hit_regions_multidom  = float(hit_info[12])
            hit_overlaps          = float(hit_info[13])
            hit_envelopes         = float(hit_info[14])
            hit_dom_n             = float(hit_info[15])
            hit_doms_within_rep_thresh = float(hit_info[16])
            hit_doms_within_inc_thresh = float(hit_info[17])
            hit_desc                   = hit_info[18]

            try:
                if hit_bitscore > high_bitscore_score[hit_seq_id]:
                    high_bitscore_score[hit_seq_id] = hit_bitscore
                    high_bitscore_line[hit_seq_id] = line
            except:
                hit_order.append(hit_seq_id)
                high_bitscore_score[hit_seq_id] = hit_bitscore
                high_bitscore_line[hit_seq_id] = line

        filtering_fields = dict()
        for hit_seq_id in hit_order:
            hit_buf.append(high_bitscore_line[hit_seq_id])
            filtering_fields[hit_seq_id] = dict()

            filter = False
            #self.log(console,"HIT_SEQ_ID: '"+hit_seq_id+"'")
            #if 'ident_thresh' in params and float(params['ident_thresh']) > float(high_bitscore_ident[hit_seq_id]):
            #    continue
            if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                filter = True
                filtering_fields[hit_seq_id]['bitscore'] = True
            #if 'overlap_fraction' in params and float(params['overlap_fraction']) > float(high_bitscore_alnlen[hit_seq_id])/float(query_len):
            #    continue
            if 'maxaccepts' in params and params['maxaccepts'] != None and hit_total == int(params['maxaccepts']):
                filter = True
                filtering_fields[hit_seq_id]['maxaccepts'] = True

            if filter:
                continue
            
            hit_total += 1
            hit_seq_ids[hit_seq_id] = True
            self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG
        

        # Measure length of hit sequences
        #
        hit_seq_len = dict()
        with open (many_forward_reads_file_path, 'r', 0) as many_forward_reads_file_handle:
            last_id = None
            last_buf = ''
            for fasta_line in many_forward_reads_file_handle.readlines():
                fasta_line = fasta_line.strip()
                if fasta_line.startswith('>'):
                    if last_id != None:
                        id_untrans = last_id
                        id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                        if id_untrans in hit_order or id_trans in hit_order:
                            hit_seq_len[last_id] = len(last_buf)
                    header = re.sub('^>', '', fasta_line)
                    last_id = re.sub('\s+.*?$', '', header)
                    last_buf = ''
                else:
                    last_buf += fasta_line
            if last_id != None:
                id_untrans = last_id
                id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                if id_untrans in hit_order or id_trans in hit_order:
                    hit_seq_len[last_id] = len(last_buf)


        # Get hit beg and end positions from Stockholm format MSA output
        #
        hit_beg = dict()
        hit_end = dict()
        with open (output_hit_MSA_file_path, 'r', 0) as output_hit_MSA_file_handle:
            for MSA_out_line in output_hit_MSA_file_handle.readlines():
                MSA_out_line = MSA_out_line.strip()
                if MSA_out_line.startswith('#=GS '):
                    hit_rec = re.sub('#=GS ', '', MSA_out_line)
                    hit_rec = re.sub('\s+.*?$', '', hit_rec)
                    hit_range = re.sub('^.*\/', '', hit_rec)
                    hit_id = re.sub('\/[^\/]+$', '', hit_rec)
                    (beg_str, end_str) = hit_range.split('-')
                    hit_beg[hit_id] = int(beg_str)
                    hit_end[hit_id] = int(end_str)


        self.log(console, 'EXTRACTING HITS FROM INPUT')
        self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


        # SequenceSet input -> SequenceSet output
        #
        if many_type_name == 'SequenceSet':
            seq_total = len(input_many_sequenceSet['sequences'])

            output_sequenceSet = dict()

            if 'sequence_set_id' in input_many_sequenceSet and input_many_sequenceSet['sequence_set_id'] != None:
                output_sequenceSet['sequence_set_id'] = input_many_sequenceSet['sequence_set_id'] + "."+search_tool_name+"_Search_filtered"
            else:
                output_sequenceSet['sequence_set_id'] = search_tool_name+"_Search_filtered"
            if 'description' in input_many_sequenceSet and input_many_sequenceSet['description'] != None:
                output_sequenceSet['description'] = input_many_sequenceSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_sequenceSet['description'] = search_tool_anme+"_Search filtered"

            self.log(console,"ADDING SEQUENCES TO SEQUENCESET")
            output_sequenceSet['sequences'] = []

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                #header_desc = seq_obj['description']
                #sequence_str = seq_obj['sequence']

                id_untrans = header_id
                id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                    #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                    accept_fids[id_untrans] = True
                    output_sequenceSet['sequences'].append(seq_obj)


        # FeatureSet input -> FeatureSet output
        #
        elif many_type_name == 'FeatureSet':
            seq_total = len(input_many_featureSet['elements'].keys())

            output_featureSet = dict()
            if 'description' in input_many_featureSet and input_many_featureSet['description'] != None:
                output_featureSet['description'] = input_many_featureSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            fId_list = input_many_featureSet['elements'].keys()
            self.log(console,"ADDING FEATURES TO FEATURESET")
            for fId in sorted(fId_list):
                for genome_ref in input_many_featureSet['elements'][fId]:
                    id_untrans = genome_ref+genome_id_feature_id_delim+fId
                    id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                    if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                        accept_fids[id_untrans] = True
                        #fId = id_untrans  # don't change fId for output FeatureSet
                        try:
                            this_genome_ref_list = output_featureSet['elements'][fId]
                        except:
                            output_featureSet['elements'][fId] = []
                            output_featureSet['element_ordering'].append(fId)
                        output_featureSet['elements'][fId].append(genome_ref)

        # Parse Genome hits into FeatureSet
        #
        elif many_type_name == 'Genome':
            seq_total = 0
            output_featureSet = dict()
#            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
#                output_featureSet['description'] = input_many_genome['scientific_name'] + " - "+search_tool_name+"_Search filtered"
#            else:
#                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()
            for fid in feature_ids:
                seq_total += 1
                id_untrans = fid
                id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                    #self.log(console, 'FOUND HIT '+fid)  # DEBUG
                    #output_featureSet['element_ordering'].append(fid)
                    accept_fids[id_untrans] = True
                    #fid = input_many_ref+genome_id_feature_id_delim+id_untrans  # don't change fId for output FeatureSet
                    output_featureSet['element_ordering'].append(fid)
                    output_featureSet['elements'][fid] = [input_many_ref]

        # Parse GenomeSet hits into FeatureSet
        #
        elif many_type_name == 'GenomeSet':
            seq_total = 0

            output_featureSet = dict()
            if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                output_featureSet['description'] = input_many_genomeSet['description'] + " - "+search_tool_name+"_Search filtered"
            else:
                output_featureSet['description'] = search_tool_name+"_Search filtered"
            output_featureSet['element_ordering'] = []
            output_featureSet['elements'] = dict()

            self.log(console,"READING HITS FOR GENOMES")  # DEBUG
            for genome_id in feature_ids_by_genome_id.keys():
                self.log(console,"READING HITS FOR GENOME "+genome_id)  # DEBUG
                genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                for feature_id in feature_ids_by_genome_id[genome_id]:
                    seq_total += 1
                    id_untrans = genome_ref+genome_id_feature_id_delim+feature_id
                    id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                    if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                        #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                        #output_featureSet['element_ordering'].append(feature['id'])
                        accept_fids[id_untrans] = True
                        #feature_id = id_untrans  # don't change fId for output FeatureSet
                        try:
                            this_genome_ref_list = output_featureSet['elements'][feature_id]
                        except:
                            output_featureSet['elements'][feature_id] = []
                            output_featureSet['element_ordering'].append(feature_id)
                        output_featureSet['elements'][feature_id].append(genome_ref)


        # load the method provenance from the context object
        #
        self.log(console,"SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
#        provenance[0]['input_ws_objects'].append(input_one_ref)
        provenance[0]['input_ws_objects'].append(input_msa_ref)
        provenance[0]['input_ws_objects'].append(input_many_ref)
        provenance[0]['service'] = 'kb_blast'
        provenance[0]['method'] = search_tool_name+'_Search'


        # Upload results
        #
        if len(invalid_msgs) == 0 and len(hit_seq_ids.keys()) > 0:
            self.log(console,"UPLOADING RESULTS")  # DEBUG

            # input many SequenceSet -> save SequenceSet
            #
            if many_type_name == 'SequenceSet':
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseSequences.SequenceSet',
                                    'data': output_sequenceSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]

            else:  # input FeatureSet, Genome, and GenomeSet -> upload FeatureSet output
                new_obj_info = ws.save_objects({
                            'workspace': params['workspace_name'],
                            'objects':[{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_featureSet,
                                    'name': params['output_filtered_name'],
                                    'meta': {},
                                    'provenance': provenance
                                }]
                        })[0]


        # build output report object
        #
        self.log(console,"BUILDING REPORT")  # DEBUG
        if len(invalid_msgs) == 0 and len(hit_order) > 0:

            # text report
            #
            report += 'sequences in search db: '+str(seq_total)+"\n"
            report += 'sequences in hit set: '+str(len(hit_order))+"\n"
            report += 'sequences in accepted hit set: '+str(hit_total)+"\n"
            report += "\n"
            for line in hit_buf:
                report += line
            self.log (console, report)


            # build html report
            if many_type_name == 'Genome':
                feature_id_to_function = GenomeToFASTA_retVal['feature_id_to_function']
                genome_ref_to_sci_name = GenomeToFASTA_retVal['genome_ref_to_sci_name']
            elif many_type_name == 'GenomeSet':
                feature_id_to_function = GenomeSetToFASTA_retVal['feature_id_to_function']
                genome_ref_to_sci_name = GenomeSetToFASTA_retVal['genome_ref_to_sci_name']
            elif many_type_name == 'FeatureSet':
                feature_id_to_function = FeatureSetToFASTA_retVal['feature_id_to_function']
                genome_ref_to_sci_name = FeatureSetToFASTA_retVal['genome_ref_to_sci_name']
                
            head_color = "#eeeeff"
            border_head_color = "#ffccff"
            accept_row_color = 'white'
            #reject_row_color = '#ffeeee'
            reject_row_color = '#eeeeee'
            reject_cell_color = '#ffcccc'
            text_fontsize = "2"
            text_color = '#606060'
            border_body_color = "#cccccc"
            bar_width = 100
            bar_height = 15
            bar_color = "lightblue"
            bar_line_color = "#cccccc"
            bar_fontsize = "1"
            bar_char = "."
            cellpadding = "3"
            cellspacing = "2"
            border = "0"

            html_report_lines = []
            html_report_lines += ['<html>']
            html_report_lines += ['<body bgcolor="white">']
            html_report_lines += ['<table cellpadding='+cellpadding+' cellspacing = '+cellspacing+' border='+border+'>']
            html_report_lines += ['<tr bgcolor="'+head_color+'">']
            html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'ALIGNMENT COVERAGE (HIT SEQ)'+'</font></td>']
            html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'GENE ID'+'</font></td>']
            html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'FUNCTION'+'</font></td>']
            html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'GENOME'+'</font></td>']
#            html_report_lines += ['<td align=center style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'IDENT'+'%</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'ALN_LEN'+'</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'E-VALUE'+'</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'BIT SCORE'+'</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'<nobr>H_BEG-H_END</nobr>'+'</font></td>']
#            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'MIS MATCH'+'</font></td>']
#            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'GAP OPEN'+'</font></td>']
            html_report_lines += ['</tr>']

            for line in hit_buf:
                line = line.strip()
                if line == '' or line.startswith('#'):
                    continue

                [hit_id, hit_accession, query_name, query_accession, e_value, bit_score, bias, e_value_best_dom, bit_score_best_dom, bias_best_dom, expected_dom_n, regions, regions_multidom, overlaps, envelopes, dom_n, doms_within_rep_thresh, doms_within_inc_thresh, hit_desc] = re.split('\s+',line)[0:19]

#                [query_id, hit_id, identity, aln_len, mismatches, gap_openings, q_beg, q_end, h_beg, h_end, e_value, bit_score] = line.split("\t")[0:12]
#                identity = str(round(float(identity), 1))
#                if identity == '100.0':  identity = '100'

                # get coords with respect to hit sequence
                h_len = hit_seq_len[hit_id]
                h_beg = hit_beg[hit_id]
                h_end = hit_end[hit_id]
                aln_len = h_end-h_beg+1
                aln_len_perc = round (100.0*float(aln_len)/float(h_len), 1)


                #if many_type_name == 'SingleEndLibrary':
                #    pass
                #elif many_type_name == 'SequenceSet':
                if many_type_name == 'SequenceSet':
                    pass
                elif many_type_name == 'Genome' or \
                        many_type_name == 'GenomeSet' or \
                        many_type_name == 'FeatureSet':

                    if many_type_name != 'Genome':
                        [genome_ref, hit_fid] = hit_id.split(genome_id_feature_id_delim)
                    else:
                        genome_ref = input_many_ref
                        hit_fid = hit_id

                    # can't just use hit_fid because may have pipes translated and can't translate back
                    fid_lookup = None
                    for fid in feature_id_to_function[genome_ref].keys():
                        id_untrans = fid
                        id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format

                        #self.log (console, "SCANNING FIDS.  HIT_FID: '"+str(hit_fid)+"' FID: '"+str(fid)+"' TRANS: '"+str(id_trans)+"'")  # DEBUG

                        if id_untrans == hit_fid or id_trans == hit_fid:
                            #self.log (console, "GOT ONE!")  # DEBUG
                            if many_type_name == 'Genome':
                                accept_id = fid
                            elif many_type_name == 'GenomeSet' or many_type_name == 'FeatureSet':
                                accept_id = genome_ref+genome_id_feature_id_delim+fid
                            if accept_id in accept_fids:
                                row_color = accept_row_color
                            else:
                                row_color = reject_row_color
                            fid_lookup = fid
                            break
                    #self.log (console, "HIT_FID: '"+str(hit_fid)+"' FID_LOOKUP: '"+str(fid_lookup)+"'")  # DEBUG
                    if fid_lookup == None:
                        raise ValueError ("unable to find fid for hit_fid: '"+str(hit_fid))
                    elif fid_lookup not in feature_id_to_function[genome_ref]:
                        raise ValueError ("unable to find function for fid: '"+str(fid_lookup))
                    fid_disp = re.sub (r"^.*\.([^\.]+)\.([^\.]+)$", r"\1.\2", fid_lookup)

                    func_disp = feature_id_to_function[genome_ref][fid_lookup]
                    genome_sci_name = genome_ref_to_sci_name[genome_ref]

                    #if 'overlap_fraction' in params and float(params['overlap_fraction']) > float(high_bitscore_alnlen[hit_seq_id])/float(query_len):

                    html_report_lines += ['<tr bgcolor="'+row_color+'">']
                    #html_report_lines += ['<tr bgcolor="'+'white'+'">']  # DEBUG
                    # add overlap bar

                    # coverage graphic (with respect to hit seq)
                    html_report_lines += ['<td valign=middle align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'">']
                    html_report_lines += ['<table style="height:'+str(bar_height)+'px; width:'+str(bar_width)+'px" border=0 cellpadding=0 cellspacing=0>']
                    full_len_pos = bar_width
                    aln_beg_pos = int (float(bar_width) * float(int(h_beg)-1)/float(int(h_len)-1))
                    aln_end_pos = int (float(bar_width) * float(int(h_end)-1)/float(int(h_len)-1))
                    cell_pix_height = str(int(round(float(bar_height)/3.0, 0)))

                    cell_color = ['','','']
                    cell_width = []
                    cell_width.append(aln_beg_pos)
                    cell_width.append(aln_end_pos-aln_beg_pos)
                    cell_width.append(bar_width-aln_end_pos)

                    for row_i in range(3):
                        html_report_lines += ['<tr style="height:'+cell_pix_height+'px">']
                        unalign_color = row_color
                        if row_i == 1:
                            unalign_color = bar_line_color
                        cell_color[0] = unalign_color
                        cell_color[1] = bar_color
                        cell_color[2] = unalign_color

                        for col_i in range(3):
                            cell_pix_width = str(cell_width[col_i])
                            cell_pix_color = cell_color[col_i]
                            html_report_lines += ['<td style="height:'+cell_pix_height+'px; width:'+cell_pix_width+'px" bgcolor="'+cell_pix_color+'"></td>']
                        html_report_lines += ['</tr>']
                    html_report_lines += ['</table>']
                    html_report_lines += ['</td>']

                    # add other cells
                    # fid
                    html_report_lines += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(fid_disp)+'</font></td>']
#                    html_report_lines += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(hit_accession)+'</font></td>']
                    # func
                    html_report_lines += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+func_disp+'</font></td>']
                    # sci name
                    html_report_lines += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+genome_sci_name+'</font></td>']
                    # ident
#                    if 'ident_thresh' in filtering_fields[hit_id]:
 #                       this_cell_color = reject_cell_color
 #                   else:
 #                       this_cell_color = row_color
 #                   html_report_lines += ['<td align=center bgcolor="'+this_cell_color+'" style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(identity)+'%</font></td>']

                    # aln len
 #                   if 'overlap_fraction' in filtering_fields[hit_id]:
 #                       this_cell_color = reject_cell_color
 #                   else:
 #                       this_cell_color = row_color
                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(aln_len)+' ('+str(aln_len_perc)+'%)</font></td>']

                    # evalue
                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'><nobr>'+str(e_value)+'</nobr></font></td>']

                    # bit score
                    if 'bitscore' in filtering_fields[hit_id]:
                        this_cell_color = reject_cell_color
                    else:
                        this_cell_color = row_color
                    html_report_lines += ['<td align=center bgcolor="'+str(this_cell_color)+'" style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'><nobr>'+str(bit_score)+'</nobr></font></td>']
                    # bias
#                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'><nobr>'+str(bias)+'</nobr><br><nobr>('+str(bias_best_dom)+')</nobr></font></td>']

                    # aln coords only for hit seq
                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'><nobr>'+str(h_beg)+'-'+str(h_end)+'</nobr></font></td>']

                    # mismatches?
#                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(mismatches)+'</font></td>']
                    # gaps?
#                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(gap_openings)+'</font></td>']

                    # regions
#                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(regions)+'</font></td>']

                    # regions_multidom
#                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(regions_multidom)+'</font></td>']

                    # overlaps
#                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(overlaps)+'</font></td>']

                    # envelopes
#                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(envelopes)+'</font></td>']

                    # expected_dom_n
#                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(expected_dom_n)+'</font></td>']

                    # doms
#                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(dom_n)+','+str(doms_within_rep_thresh)+','+str(doms_within_inc_thresh)+'</font></td>']

                    # hit desc
#                    html_report_lines += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(hit_desc)+'</font></td>']

                    html_report_lines += ['</tr>']

            html_report_lines += ['</table>']
            html_report_lines += ['</body>']
            html_report_lines += ['</html>']


            # write html to file and upload
            #
            html_report_str = "\n".join(html_report_lines)
            html_output_dir = os.path.join(self.output_dir,'html_output')
            if not os.path.exists(html_output_dir):
                os.makedirs(html_output_dir)
            html_file = search_tool_name+'_Search.html'
            html_path = os.path.join (html_output_dir, html_file)
            with open (html_path, 'w', 0) as html_handle:
                html_handle.write(html_report_str)

            dfu = DFUClient(self.callbackURL)
            try:
                #HTML_upload_ret = dfu.file_to_shock({'file_path': html_path,
                HTML_upload_ret = dfu.file_to_shock({'file_path': html_output_dir,
                                                     'make_handle': 0,
                                                     'pack': 'zip'})
            except:
                raise ValueError ('Logging exception loading HTML file to shock')
            try:
                TAB_upload_ret = dfu.file_to_shock({'file_path': output_hit_TAB_file_path,
                                                    'make_handle': 0})
                                                    #'pack': 'zip'})
            except:
                raise ValueError ('Logging exception loading TAB output to shock')
            try:
                MSA_upload_ret = dfu.file_to_shock({'file_path': output_hit_MSA_file_path,
                                                    'make_handle': 0})
                                                    #'pack': 'zip'})
            except:
                raise ValueError ('Logging exception loading MSA output to shock')


            # create report object
            reportName = 'hmmer_report_'+str(uuid.uuid4())
            reportObj = {'objects_created': [],
                         #'text_message': '',  # or is it 'message'?
                         'message': '',  # or is it 'text_message'?
                         'direct_html': None,
                         'direct_html_link_index': None,
                         'file_links': [],
                         'html_links': [],
                         'workspace_name': params['workspace_name'],
                         'report_object_name': reportName
                         }
            #html_buf_lim = 16000  # really 16KB, but whatever
            #if len(html_report_str) <= html_buf_lim:
            #    reportObj['direct_html'] = html_report_str
            #else:
            reportObj['direct_html_link_index'] = 0
            reportObj['html_links'] = [{'shock_id': HTML_upload_ret['shock_id'],
                                        'name': html_file,
                                        'label': search_tool_name+' HTML Report'}
                                       ]

            reportObj['file_links'] = [{'shock_id': TAB_upload_ret['shock_id'],
                                        'name': search_tool_name+'_Search.TAB',
                                        'label': search_tool_name+' hits TABLE'},

                                       {'shock_id': MSA_upload_ret['shock_id'],
                                        'name': search_tool_name+'_Search.MSA',
                                        'label': search_tool_name+' hits MSA'},

                                       ]
#            if extra_output:
#                extension = 'txt'
#                if params['output_extra_format'] == '5':
#                    extension = 'xml'
#                elif params['output_extra_format'] == '8':
#                    extension = 'asn1txt'
#                elif params['output_extra_format'] == '9':
#                    extension = 'asn1bin'
#                elif params['output_extra_format'] == '10':
#                    extension = 'csv'
#                elif params['output_extra_format'] == '11':
#                    extension = 'asn1arc'
#                reportObj['file_links'].append({'shock_id': extra_upload_ret['shock_id'],
#                                                'name': search_tool_name+'_Search-m'+str(params['output_extra_format'])+'.'+extension,
#                                                'label': search_tool_name+' Results: m'+str(params['output_extra_format'])})
                            
            reportObj['objects_created'].append({'ref':str(params['workspace_name'])+'/'+params['output_filtered_name'],'description':search_tool_name+' hits'})
            #reportObj['message'] = report


            # save report object
            #
            SERVICE_VER = 'release'
            reportClient = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
            #report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})
            report_info = reportClient.create_extended_report(reportObj)

        else:
            if len(hit_order) == 0:  # no hits
                report += "No hits were found\n"
            else:  # data validation error
                report += "FAILURE\n\n"+"\n".join(invalid_msgs)+"\n"

            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            reportName = 'hmmer_report_'+str(uuid.uuid4())
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
                            'provenance':provenance
                            }
                        ]
                    })[0]
            report_info = dict()
            report_info['name'] = report_obj_info[1]
            report_info['ref'] = str(report_obj_info[6])+'/'+str(report_obj_info[0])+'/'+str(report_obj_info[4])

        self.log(console,"BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = { 'report_name': report_info['name'],
                      'report_ref': report_info['ref']
                      }
        self.log(console,search_tool_name+"_Search DONE")
        #END HMMER_MSA_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method HMMER_MSA_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def HMMER_Local_MSA_Group_Search(self, ctx, params):
        """
        Methods for HMMER search of a Local MSA Group (found automatically within workspace) against many sequences 
        **
        **    overloading as follows:
        **        input_many_ref: SingleEndLibrary, FeatureSet, Genome, GenomeSet
        **        output_name: SingleEndLibrary (if input_many is SELib), (else) FeatureSet
        :param params: instance of type "HMMER_Local_MSA_Group_Params" (HMMER
           Local MSA Group Input Params) -> structure: parameter
           "workspace_name" of type "workspace_name" (** The workspace object
           refs are of form: ** **    objects = ws.get_objects([{'ref':
           params['workspace_id']+'/'+params['obj_name']}]) ** ** "ref" means
           the entire name combining the workspace id and the object name **
           "id" is a numerical identifier of the workspace or object, and
           should just be used for workspace ** "name" is a string identifier
           of a workspace or object.  This is received from Narrative.),
           parameter "input_msa_refs" of type "data_obj_ref", parameter
           "input_many_ref" of type "data_obj_ref", parameter
           "output_filtered_name" of type "data_obj_name", parameter
           "coalesce_output" of type "bool", parameter "e_value" of Double,
           parameter "bitscore" of Double, parameter "maxaccepts" of Double,
           parameter "heatmap" of type "bool", parameter "vertical" of type
           "bool", parameter "show_blanks" of type "bool"
        :returns: instance of type "HMMER_Output" (HMMER Output) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN HMMER_Local_MSA_Group_Search
        console = []
        invalid_msgs = []
        search_tool_name = 'HMMER_Local_MSA_Group_prot'
        self.log(console,'Running '+search_tool_name+'_Search with params=')
        self.log(console, "\n"+pformat(params))
        report = ''
#        report = 'Running '+search_tool_name+'_Search with params='
#        report += "\n"+pformat(params)
        #appropriate_sequence_found_in_one_input = False
        appropriate_sequence_found_in_MSA_input = False
        appropriate_sequence_found_in_many_input = False
        genome_id_feature_id_delim = '.f:'


        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
#        if 'input_one_ref' not in params:
#            raise ValueError('input_one_ref parameter is required')
#        if 'input_msa_refs' not in params or len(params['input_msa_refs']) == 0:
#            raise ValueError('input_msa_refs parameter is required if selecting local MSAs')
        if 'input_many_ref' not in params:
            raise ValueError('input_many_ref parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')
        if 'coalesce_output' not in params:
            raise ValueError('coalesce_output parameter is required')


        # set local names and ids
#        input_one_ref = params['input_one_ref']
        #input_msa_ref = params['input_msa_ref']
        input_many_ref = params['input_many_ref']
        ws_id = input_many_ref.split('/')[0]
        

        #### Get the input_many object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            #objects = ws.get_objects([{'ref': input_many_ref}])
            objects = ws.get_objects2({'objects':[{'ref': input_many_ref}]})['data']
            input_many_data = objects[0]['data']
            info = objects[0]['info']
            input_many_name = str(info[1])
            many_type_name = info[2].split('.')[1].split('-')[0]

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()


        # Handle overloading (input_many can be SequenceSet, FeatureSet, Genome, or GenomeSet)
        #
        if many_type_name == 'SequenceSet':
            try:
                input_many_sequenceSet = input_many_data
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to get SequenceSet: ' + str(e))

            header_id = input_many_sequenceSet['sequences'][0]['sequence_id']
            many_forward_reads_file_path = os.path.join(self.output_dir, header_id+'.fasta')
            many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
            self.log(console, 'writing reads file: '+str(many_forward_reads_file_path))

            for seq_obj in input_many_sequenceSet['sequences']:
                header_id = seq_obj['sequence_id']
                sequence_str = seq_obj['sequence']

                PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
                #DNA_pattern = re.compile("^[acgtuACGTUnryNRY ]+$")
                if not PROT_pattern.match(sequence_str):
                    self.log(invalid_msgs,"BAD record for sequence_id: "+header_id+"\n"+sequence_str+"\n")
                    continue
                appropriate_sequence_found_in_many_input = True
                many_forward_reads_file_handle.write('>'+header_id+"\n")
                many_forward_reads_file_handle.write(sequence_str+"\n")
            many_forward_reads_file_handle.close();
            self.log(console, 'done')


        # FeatureSet
        #
        if many_type_name == 'FeatureSet':
            # retrieve sequences for features
            input_many_featureSet = input_many_data
            many_forward_reads_file_dir = self.output_dir
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            FeatureSetToFASTA_params = {
                'featureSet_ref':      input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            #SERVICE_VER = 'release'
            SERVICE_VER = 'dev'
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
            FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA (FeatureSetToFASTA_params)
            many_forward_reads_file_path = FeatureSetToFASTA_retVal['fasta_file_path']
            feature_ids_by_genome_ref = FeatureSetToFASTA_retVal['feature_ids_by_genome_ref']
            if len(feature_ids_by_genome_ref.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Genome
        #
        elif many_type_name == 'Genome':
            many_forward_reads_file_dir = self.output_dir
            many_forward_reads_file = input_many_name+".fasta"

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeToFASTA_params = {
                'genome_ref':          input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            #SERVICE_VER = 'release'
            SERVICE_VER = 'dev'
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
            GenomeToFASTA_retVal = DOTFU.GenomeToFASTA (GenomeToFASTA_params)
            many_forward_reads_file_path = GenomeToFASTA_retVal['fasta_file_path']
            feature_ids = GenomeToFASTA_retVal['feature_ids']
            if len(feature_ids) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "Genome2Fasta() took "+str(end_time-beg_time)+" secs")
            

        # GenomeSet
        #
        elif many_type_name == 'GenomeSet':
            input_many_genomeSet = input_many_data
            many_forward_reads_file_dir = self.output_dir
            many_forward_reads_file = input_many_name+".fasta"

            genome_refs = []
            for genome_id in input_many_genomeSet['elements']:
                genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                if genome_ref not in genome_refs:
                    genome_refs.append(genome_ref)

            # DEBUG
            #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            GenomeSetToFASTA_params = {
                'genomeSet_ref':       input_many_ref,
                'file':                many_forward_reads_file,
                'dir':                 many_forward_reads_file_dir,
                'console':             console,
                'invalid_msgs':        invalid_msgs,
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }

            #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
            #SERVICE_VER = 'release'
            SERVICE_VER = 'dev'
            DOTFU = KBaseDataObjectToFileUtils (url=self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
            GenomeSetToFASTA_retVal = DOTFU.GenomeSetToFASTA (GenomeSetToFASTA_params)
            many_forward_reads_file_path = GenomeSetToFASTA_retVal['fasta_file_path_list'][0]
            feature_ids_by_genome_id = GenomeSetToFASTA_retVal['feature_ids_by_genome_id']
            if len(feature_ids_by_genome_id.keys()) > 0:
                appropriate_sequence_found_in_many_input = True

            # DEBUG
            #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
            #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")


        # Missing proper input_many_type
        #
        else:
            raise ValueError('Cannot yet handle input_many type of: '+many_type_name)            


        #### Get the input_msa_refs
        ##
        if 'input_msa_refs' in params and len(params['input_msa_refs']) != 0:
            input_msa_refs = params['input_msa_refs']
        else:
            input_msa_refs = []
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            try:
                msa_obj_info_list = ws.list_objects({'ids':[ws_id],'type':"KBaseTrees.MSA"})
            except Exception as e:
                raise ValueError ("Unable to list MSA objects from workspace: "+str(params['workspace_name'])+" "+str(e))

            for info in msa_obj_info_list:
                [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
                
                input_msa_ref = str(info[WSID_I])+'/'+str(info[OBJID_I])+'/'+str(info[VERSION_I])
                input_msa_refs.append(input_msa_ref)


        #### write the MSAs to file and collect names and desc
        ##
        input_msa_names = []
        input_msa_descs = []
        for input_msa_ref in input_msa_refs:
            try:
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
                #objects = ws.get_objects([{'ref': input_msa_ref}])
                objects = ws.get_objects2({'objects':[{'ref': input_msa_ref}]})['data']
                input_msa_data = objects[0]['data']
                info = objects[0]['info']
                input_msa_name = str(info[1])
                msa_type_name = info[2].split('.')[1].split('-')[0]

            except Exception as e:
                raise ValueError('Unable to fetch '+input_msa_name+' object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()

            # set hmmer_dir
            hmmer_dir = os.path.join(self.output_dir, input_msa_name)
            if not os.path.exists(hmmer_dir):
                os.makedirs(hmmer_dir)

            if msa_type_name == 'MSA':
                self.log (console, "\n\nPROCESSING MSA "+input_msa_name+"\n")  # DEBUG

                input_msa_names.append(input_msa_name)

                MSA_in = input_msa_data
                if 'description' in MSA_in and MSA_in['description'] != None and MSA_in['description'] != '':
                    input_msa_descs.append(MSA_in['description'])
                else:
                    input_msa_descs.append(input_msa_name)

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
                input_MSA_file_path = os.path.join(hmmer_dir, input_msa_name+".clustal")
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
                conservation_symbol = ''
                for i in range(row_len):
                    first_seen_char = MSA_in['alignment'][row_id_0][i]
                    symbol = '*'
                    for row_id in row_order:
                        if MSA_in['alignment'][row_id][i] == '-' or MSA_in['alignment'][row_id][i] != first_seen_char:
                            symbol = ' '
                            break
                    conservation_symbol += symbol

                # break up MSA into 60 char chunks
                records = []
                chunk_len = 60
                whole_chunks = int(math.floor(row_len/chunk_len))
                if whole_chunks > 0:
                    for j in range(whole_chunks):
                        records.append('')
                        for row_id in row_order:
                            padding = ''
                            if longest_row_id_len-len(row_id) > 0:
                                for i in range(0,longest_row_id_len-len(row_id)):
                                    padding += ' '
                            records.append(row_id + padding + " " +
                                       MSA_in['alignment'][row_id][j*chunk_len:(j+1)*chunk_len])
                        records.append(''.join([' ' for s in range(longest_row_id_len)]) + " " +
                                   conservation_symbol[j*chunk_len:(j+1)*chunk_len])

                # add final rows
                if (row_len % chunk_len) != 0:
                    j=whole_chunks
                    records.append('')
                    for row_id in row_order:
                        padding = ''
                        if longest_row_id_len-len(row_id) > 0:
                            for i in range(0,longest_row_id_len-len(row_id)):
                                padding += ' '
                        records.append(row_id + padding + " " +
                                   MSA_in['alignment'][row_id][j*chunk_len:row_len])
                    records.append(''.join([' ' for s in range(longest_row_id_len)]) + " " +
                               conservation_symbol[j*chunk_len:row_len])
            
                # write that sucker
                with open(input_MSA_file_path,'w',0) as input_MSA_file_handle:
                    input_MSA_file_handle.write(header+"\n")
                    input_MSA_file_handle.write("\n".join(records)+"\n")

                # DEBUG
                #report += "MSA:\n"
                #report += header+"\n"
                #report += "\n".join(records)+"\n"
                #self.log(console,report)


                # Determine whether nuc or protein sequences
                #
                self.log (console, "CHECKING MSA for PROTEIN seqs...")  # DEBUG
                PROT_MSA_pattern = re.compile("^[\.\-_acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
                #NUC_MSA_pattern = re.compile("^[\.\-_ACGTUXNRYSWKMBDHVacgtuxnryswkmbdhv \t\n]+$")
                appropriate_sequence_found_in_MSA_input = True
                for row_id in row_order:
                    #self.log(console, row_id+": '"+MSA_in['alignment'][row_id]+"'")    # DEBUG
                    if not PROT_MSA_pattern.match(MSA_in['alignment'][row_id]):
                        self.log(invalid_msgs,"BAD record for MSA row_id: "+row_id+"\n"+MSA_in['alignment'][row_id]+"\n")
                        appropriate_sequence_found_in_MSA_input = False
                        break

                if not appropriate_sequence_found_in_MSA_input:
                    self.log(invalid_msgs,"no protein sequences found in '"+input_msa_name+"'")


            # Missing proper input_type
            #
            else:
                raise ValueError('Cannot yet handle input_msa type of: '+msa_type_name)


        # check for failed input file creation
        #
#        if not appropriate_sequence_found_in_one_input:
#            self.log(invalid_msgs,"no protein sequences found in '"+input_one_name+"'")
        if not appropriate_sequence_found_in_many_input:
            self.log(invalid_msgs,"no protein sequences found in '"+input_many_name+"'")


        # input data failed validation.  Need to return
        #
        if len(invalid_msgs) > 0:

            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
#            provenance[0]['input_ws_objects'].append(input_one_ref)
            provenance[0]['input_ws_objects'].append(input_many_ref)
            for input_msa_ref in input_msa_refs:
                provenance[0]['input_ws_objects'].append(input_msa_ref)
            provenance[0]['service'] = 'kb_hmmer'
            provenance[0]['method'] = search_tool_name+'_Search'


            # build output report object
            #
            self.log(console,"BUILDING REPORT")  # DEBUG
            report += "FAILURE:\n\n"+"\n".join(invalid_msgs)+"\n"
            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            reportName = 'hmmer_report_'+str(uuid.uuid4())
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            report_obj_info = ws.save_objects({
                    #'id':info[6],
                    'workspace':params['workspace_name'],
                    'objects':[
                        {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance  # DEBUG
                        }
                        ]
                    })[0]

            self.log(console,"BUILDING RETURN OBJECT")
            returnVal = { 'report_name': reportName,
                      'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      }
            self.log(console,search_tool_name+"_Search DONE")
            return [returnVal]

        
        #### iterate through MSAs and scan input_many DBs
        ##
        total_hit_cnt = []
        output_hit_TAB_file_paths = []
        output_hit_MSA_file_paths = []
        output_filtered_fasta_file_paths = []
        output_hits_flags = []
        objects_created_refs = []
        coalesced_sequenceObjs = []
        coalesce_featureIds_element_ordering = []
        coalesce_featureIds_genome_ordering = []
        html_report_chunks = []
        hit_cnt_by_genome_and_model = dict()


        for i,input_msa_ref in enumerate(input_msa_refs):

            ### set paths
            #
            input_msa_name = input_msa_names[i]
            hmmer_dir = os.path.join(self.output_dir, input_msa_name)  # this must match above
            input_MSA_file_path = os.path.join(hmmer_dir, input_msa_name+".clustal")

            #output_aln_file_path = os.path.join(hmmer_dir, input_msa_name+'.alnout.txt');
            #output_extra_file_path = os.path.join(hmmer_dir, input_msa_name+'.alnout_extra.txt');
            #output_filtered_fasta_file_path = os.path.join(hmmer_dir, input_msa_name+'.output_filtered.faa');


            ### Build HMM from MSA
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
            hmmer_build_cmd.append('CLUSTAL')
            hmmer_build_cmd.append(HMM_file_path)
            hmmer_build_cmd.append(input_MSA_file_path)

            # Run HMMER_BUILD, capture output as it happens
            #
            self.log(console, 'RUNNING HMMER_BUILD:')
            self.log(console, '    '+' '.join(hmmer_build_cmd))
            #report += "\n"+'running HMMER_BUILD:'+"\n"
            #report += '    '+' '.join(hmmer_build_cmd)+"\n"

            p = subprocess.Popen(hmmer_build_cmd, \
                                     cwd = self.output_dir, \
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
            elif not os.path.getsize(HMM_file_path) > 0:
                raise ValueError("HMMER_BUILD created empty HMM file '"+HMM_file_path+"'")


            ### Construct the HMMER_SEARCH command
            #
            # SYNTAX (from http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf)
            #
            # hmmsearch --tblout <TAB_out> -A <MSA_out> --noali --notextw -E <e_value> -T <bit_score> <hmmfile> <seqdb>
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

            output_hit_TAB_file_path = os.path.join(hmmer_dir, input_msa_name+'.hitout.txt');
            output_hit_MSA_file_path = os.path.join(hmmer_dir, input_msa_name+'.msaout.txt');
            output_filtered_fasta_file_path = os.path.join(hmmer_dir, input_msa_name+'.output_filtered.fasta');
            output_hit_TAB_file_paths.append (output_hit_TAB_file_path)
            output_hit_MSA_file_paths.append (output_hit_MSA_file_path)
            output_filtered_fasta_file_paths.append(output_filtered_fasta_file_path)


            # this is command for basic search mode
            hmmer_search_cmd.append('--tblout')
            hmmer_search_cmd.append(output_hit_TAB_file_path)
            hmmer_search_cmd.append('-A')
            hmmer_search_cmd.append(output_hit_MSA_file_path)
            hmmer_search_cmd.append('--noali')
            hmmer_search_cmd.append('--notextw')
            hmmer_search_cmd.append('-E')  # can't use -T with -E, so we'll use -E
            hmmer_search_cmd.append(str(params['e_value']))
            hmmer_search_cmd.append(HMM_file_path)
            hmmer_search_cmd.append(many_forward_reads_file_path)

            # options
            #if 'maxaccepts' in params:
            #    if params['maxaccepts']:
            #        hmmer_search_cmd.append('-max_target_seqs')
            #        hmmer_search_cmd.append(str(params['maxaccepts']))
            
            # Run HMMER, capture output as it happens
            #
            self.log(console, 'RUNNING HMMER_SEARCH:')
            self.log(console, '    '+' '.join(hmmer_search_cmd))
            #report += "\n"+'running HMMER_SEARCH:'+"\n"
            #report += '    '+' '.join(hmmer_search_cmd)+"\n"

            p = subprocess.Popen(hmmer_search_cmd, \
                                     cwd = self.output_dir, \
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
            elif not os.path.getsize(output_hit_TAB_file_path) > 0:
                raise ValueError("HMMER_SEARCH created empty TAB file '"+output_hit_TAB_file_path+"'")
            if not os.path.isfile(output_hit_MSA_file_path):
                raise ValueError("HMMER_SEARCH failed to create MSA file '"+output_hit_MSA_file_path+"'")
            elif not os.path.getsize(output_hit_MSA_file_path) > 0:
                #raise ValueError("HMMER_SEARCH created empty MSA file '"+output_hit_MSA_file_path+"'")
                self.log(console,"HMMER_SEARCH created empty MSA file '"+output_hit_MSA_file_path+"'")
                objects_created_refs.append(None)
                total_hit_cnt.append(0)
                continue


            # DEBUG
            #self.log(console, "DEBUG: output_hit_TAB_file_path: '"+str(output_hit_TAB_file_path))
            #self.log(console, "DEBUG: output_hit_MSA_file_path: '"+str(output_hit_MSA_file_path))
            #report = "TAB:\n\n"
            #with open (output_hit_TAB_file_path, 'r') as output_handle:
            #    for line in output_handle:
            #        report += line+"\n"
            #report += "\n\nMSA:\n\n"
            #with open (output_hit_MSA_file_path, 'r') as output_handle:
            #    for line in output_handle:
            #        report += line+"\n"
            #self.log(console, report)


            ### Parse the HMMER tabular output and store ids to filter many set to make filtered object to save back to KBase
            #
            self.log(console, 'PARSING HMMER SEARCH OUTPUT')
            if not os.path.isfile(output_hit_TAB_file_path):
                raise ValueError("failed to create HMMER output: "+output_hit_TAB_file_path)
            elif not os.path.getsize(output_hit_TAB_file_path) > 0:
                raise ValueError("created empty file for HMMER output: "+output_hit_TAB_file_path)
            hit_seq_ids = dict()
            accept_fids = dict()
            output_hit_TAB_file_handle = open (output_hit_TAB_file_path, "r", 0)
            output_aln_buf = output_hit_TAB_file_handle.readlines()
            output_hit_TAB_file_handle.close()
            hit_total = 0
            high_bitscore_line = dict()
            high_bitscore_score = dict()
            high_bitscore_ident = dict()
            high_bitscore_alnlen = dict()
            hit_order = []
            hit_buf = []
            hit_accept_something = False
            #header_done = False
            for line in output_aln_buf:
                if line.startswith('#'):
                    #if not header_done:
                    #    hit_buf.append(line)
                    continue
                #header_done = True
                #self.log(console,'HIT LINE: '+line)  # DEBUG
                hit_info = re.split ('\s+', line)
                hit_seq_id            = hit_info[0]
                hit_accession         = hit_info[1]
                query_name            = hit_info[2]
                query_accession       = hit_info[3]
                hit_e_value           = float(hit_info[4])
                hit_bitscore          = float(hit_info[5])
                hit_bias              = float(hit_info[6])
                hit_e_value_best_dom  = float(hit_info[7])
                hit_bitscore_best_dom = float(hit_info[8])
                hit_bias_best_dom     = float(hit_info[9])
                hit_expected_dom_n    = float(hit_info[10])
                hit_regions           = float(hit_info[11])
                hit_regions_multidom  = float(hit_info[12])
                hit_overlaps          = float(hit_info[13])
                hit_envelopes         = float(hit_info[14])
                hit_dom_n             = float(hit_info[15])
                hit_doms_within_rep_thresh = float(hit_info[16])
                hit_doms_within_inc_thresh = float(hit_info[17])
                hit_desc                   = hit_info[18]

                try:
                    if hit_bitscore > high_bitscore_score[hit_seq_id]:
                        high_bitscore_score[hit_seq_id] = hit_bitscore
                        high_bitscore_line[hit_seq_id] = line
                except:
                    hit_order.append(hit_seq_id)
                    high_bitscore_score[hit_seq_id] = hit_bitscore
                    high_bitscore_line[hit_seq_id] = line

            filtering_fields = dict()
            total_hit_cnt.append (len(hit_order))

            for hit_seq_id in hit_order:
                hit_buf.append(high_bitscore_line[hit_seq_id])
                filtering_fields[hit_seq_id] = dict()

                filter = False
                #self.log(console,"HIT_SEQ_ID: '"+hit_seq_id+"'")
                #if 'ident_thresh' in params and float(params['ident_thresh']) > float(high_bitscore_ident[hit_seq_id]):
                #    continue
                if 'bitscore' in params and float(params['bitscore']) > float(high_bitscore_score[hit_seq_id]):
                    filter = True
                    filtering_fields[hit_seq_id]['bitscore'] = True
                #if 'overlap_fraction' in params and float(params['overlap_fraction']) > float(high_bitscore_alnlen[hit_seq_id])/float(query_len):
                #    continue
                if 'maxaccepts' in params and params['maxaccepts'] != None and hit_total == int(params['maxaccepts']):
                    filter = True
                    filtering_fields[hit_seq_id]['maxaccepts'] = True

                if filter:
                    continue
            
                hit_accept_something = True
                hit_total += 1
                hit_seq_ids[hit_seq_id] = True
                self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG

                # capture accepted hit count by genome_ref and model
                genome_ref = hit_seq_id.split(genome_id_feature_id_delim)[0]
                self.log(console, "DEBUG: genome_ref: '"+str(genome_ref)+"'")
                self.log(console, "DEBUG: input_msa_name: '"+str(input_msa_name)+"'")
                if genome_ref not in hit_cnt_by_genome_and_model:
                    hit_cnt_by_genome_and_model[genome_ref] = dict()
                if input_msa_name not in hit_cnt_by_genome_and_model[genome_ref]:
                    hit_cnt_by_genome_and_model[genome_ref][input_msa_name] = 0
                hit_cnt_by_genome_and_model[genome_ref][input_msa_name] += 1


            # Measure length of hit sequences
            #
            hit_seq_len = dict()
            with open (many_forward_reads_file_path, 'r', 0) as many_forward_reads_file_handle:
                last_id = None
                last_buf = ''
                for fasta_line in many_forward_reads_file_handle.readlines():
                    fasta_line = fasta_line.strip()
                    if fasta_line.startswith('>'):
                        if last_id != None:
                            id_untrans = last_id
                            id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                            if id_untrans in hit_order or id_trans in hit_order:
                                hit_seq_len[last_id] = len(last_buf)
                        header = re.sub('^>', '', fasta_line)
                        last_id = re.sub('\s+.*?$', '', header)
                        last_buf = ''
                    else:
                        last_buf += fasta_line
                if last_id != None:
                    id_untrans = last_id
                    id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                    if id_untrans in hit_order or id_trans in hit_order:
                        hit_seq_len[last_id] = len(last_buf)


            # Get hit beg and end positions from Stockholm format MSA output
            #
            hit_beg = dict()
            hit_end = dict()
            with open (output_hit_MSA_file_path, 'r', 0) as output_hit_MSA_file_handle:
                for MSA_out_line in output_hit_MSA_file_handle.readlines():
                    MSA_out_line = MSA_out_line.strip()
                    if MSA_out_line.startswith('#=GS '):
                        hit_rec = re.sub('#=GS ', '', MSA_out_line)
                        hit_rec = re.sub('\s+.*?$', '', hit_rec)
                        hit_range = re.sub('^.*\/', '', hit_rec)
                        hit_id = re.sub('\/[^\/]+$', '', hit_rec)
                        (beg_str, end_str) = hit_range.split('-')
                        hit_beg[hit_id] = int(beg_str)
                        hit_end[hit_id] = int(end_str)

            #
            ### Create output objects
            #

            self.log(console, 'EXTRACTING HITS FROM INPUT')
            self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG


            # SequenceSet input -> SequenceSet output
            #
            if many_type_name == 'SequenceSet':
                seq_total = len(input_many_sequenceSet['sequences'])

                output_sequenceSet = dict()

                if 'sequence_set_id' in input_many_sequenceSet and input_many_sequenceSet['sequence_set_id'] != None:
                    output_sequenceSet['sequence_set_id'] = input_many_sequenceSet['sequence_set_id'] + "."+search_tool_name+"_Search_filtered"
                else:
                    output_sequenceSet['sequence_set_id'] = search_tool_name+"_Search_filtered"
                if 'description' in input_many_sequenceSet and input_many_sequenceSet['description'] != None:
                    output_sequenceSet['description'] = input_many_sequenceSet['description'] + " - "+search_tool_name+"_Search filtered"
                else:
                    output_sequenceSet['description'] = search_tool_anme+"_Search filtered"

                self.log(console,"ADDING SEQUENCES TO SEQUENCESET")
                output_sequenceSet['sequences'] = []

                for seq_obj in input_many_sequenceSet['sequences']:
                    header_id = seq_obj['sequence_id']
                    #header_desc = seq_obj['description']
                    #sequence_str = seq_obj['sequence']

                    id_untrans = header_id
                    id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                    if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                        #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                        accept_fids[id_untrans] = True
                        output_sequenceSet['sequences'].append(seq_obj)


            # FeatureSet input -> FeatureSet output
            #
            elif many_type_name == 'FeatureSet':
                seq_total = len(input_many_featureSet['elements'].keys())

                output_featureSet = dict()
                if 'description' in input_many_featureSet and input_many_featureSet['description'] != None:
                    output_featureSet['description'] = input_many_featureSet['description'] + " - "+search_tool_name+"_Search filtered"
                else:
                    output_featureSet['description'] = search_tool_name+"_Search filtered"
                output_featureSet['element_ordering'] = []
                output_featureSet['elements'] = dict()

                fId_list = input_many_featureSet['elements'].keys()
                self.log(console,"ADDING FEATURES TO FEATURESET")
                for fId in sorted(fId_list):
                    for genome_ref in input_many_featureSet['elements'][fId]:
                        id_untrans = genome_ref+genome_id_feature_id_delim+fId
                        id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                        if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                            #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                            accept_fids[id_untrans] = True
                            #fId = id_untrans  # don't change fId for output FeatureSet
                            try:
                                this_genome_ref_list = output_featureSet['elements'][fId]
                            except:
                                output_featureSet['elements'][fId] = []
                                output_featureSet['element_ordering'].append(fId)
                            output_featureSet['elements'][fId].append(genome_ref)

            # Parse Genome hits into FeatureSet
            #
            elif many_type_name == 'Genome':
                seq_total = 0
                output_featureSet = dict()
    #            if 'scientific_name' in input_many_genome and input_many_genome['scientific_name'] != None:
    #                output_featureSet['description'] = input_many_genome['scientific_name'] + " - "+search_tool_name+"_Search filtered"
    #            else:
    #                output_featureSet['description'] = search_tool_name+"_Search filtered"
                output_featureSet['description'] = search_tool_name+"_Search filtered"
                output_featureSet['element_ordering'] = []
                output_featureSet['elements'] = dict()
                for fid in feature_ids:
                    seq_total += 1
                    id_untrans = fid
                    id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                    if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                        #self.log(console, 'FOUND HIT '+fid)  # DEBUG
                        #output_featureSet['element_ordering'].append(fid)
                        accept_fids[id_untrans] = True
                        #fid = input_many_ref+genome_id_feature_id_delim+id_untrans  # don't change fId for output FeatureSet
                        output_featureSet['element_ordering'].append(fid)
                        output_featureSet['elements'][fid] = [input_many_ref]

            # Parse GenomeSet hits into FeatureSet
            #
            elif many_type_name == 'GenomeSet':
                seq_total = 0

                output_featureSet = dict()
                if 'description' in input_many_genomeSet and input_many_genomeSet['description'] != None:
                    output_featureSet['description'] = input_many_genomeSet['description'] + " - "+search_tool_name+"_Search filtered"
                else:
                    output_featureSet['description'] = search_tool_name+"_Search filtered"
                output_featureSet['element_ordering'] = []
                output_featureSet['elements'] = dict()

                self.log(console,"READING HITS FOR GENOMES")  # DEBUG
                for genome_id in feature_ids_by_genome_id.keys():
                    self.log(console,"READING HITS FOR GENOME "+genome_id)  # DEBUG
                    genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                    for feature_id in feature_ids_by_genome_id[genome_id]:
                        seq_total += 1
                        id_untrans = genome_ref+genome_id_feature_id_delim+feature_id
                        id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format
                        if id_trans in hit_seq_ids or id_untrans in hit_seq_ids:
                            #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                            #output_featureSet['element_ordering'].append(feature['id'])
                            accept_fids[id_untrans] = True
                            #feature_id = id_untrans  # don't change fId for output FeatureSet
                            try:
                                this_genome_ref_list = output_featureSet['elements'][feature_id]
                            except:
                                output_featureSet['elements'][feature_id] = []
                                output_featureSet['element_ordering'].append(feature_id)
                            output_featureSet['elements'][feature_id].append(genome_ref)


            # load the method provenance from the context object
            #
            self.log(console,"SETTING PROVENANCE")  # DEBUG
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            # add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects'] = []
    #        provenance[0]['input_ws_objects'].append(input_one_ref)
            provenance[0]['input_ws_objects'].append(input_msa_ref)
            provenance[0]['input_ws_objects'].append(input_many_ref)
            provenance[0]['service'] = 'kb_blast'
            provenance[0]['method'] = search_tool_name+'_Search'


            ### Create output object
            #
            if 'coalesce_output' in params and int(params['coalesce_output']) == 1:
                if len(invalid_msgs) == 0:
                    if len(hit_seq_ids.keys()) == 0:   # Note, this is after filtering, so there may be more unfiltered hits
                        self.log(console,"No Object to Upload for MSA "+input_msa_name)  # DEBUG
                        objects_created_refs.append(None)
                        continue

                    # accumulate hits into coalesce object
                    #
                    if many_type_name == 'SequenceSet':  # input many SequenceSet -> save SequenceSet
                        for seq_obj in output_sequenceSet['sequences']:
                            coalesced_sequenceObjs.append(seq_obj)

                    else:  # input FeatureSet, Genome, and GenomeSet -> upload FeatureSet output

                        for fId in output_featureSet['element_ordering']:
                            coalesce_featureIds_element_ordering.append(fId)
                            #coalesce_featureIds_genome_ordering.append(output_featureSet['elements'][fId][0])
                            for this_genome_ref in output_featureSet['elements'][fId]:
                                coalesce_featureIds_genome_ordering.append(this_genome_ref)

            else:  # keep output separate  Upload results if coalesce_output is 0
                output_name = input_msa_name+'-'+params['output_filtered_name']

                if len(invalid_msgs) == 0:
                    if len(hit_seq_ids.keys()) == 0:   # Note, this is after filtering, so there may be more unfiltered hits
                        self.log(console,"No Object to Upload for MSA "+input_msa_name)  # DEBUG
                        objects_created_refs.append(None)
                        continue

                    self.log(console,"Uploading results Object MSA "+input_msa_name)  # DEBUG

                    # input many SequenceSet -> save SequenceSet
                    #
                    if many_type_name == 'SequenceSet':
                        new_obj_info = ws.save_objects({
                                    'workspace': params['workspace_name'],
                                    'objects':[{
                                            'type': 'KBaseSequences.SequenceSet',
                                            'data': output_sequenceSet,
                                            'name': output_name,
                                            'meta': {},
                                            'provenance': provenance
                                        }]
                                })[0]

                    else:  # input FeatureSet, Genome, and GenomeSet -> upload FeatureSet output
                        new_obj_info = ws.save_objects({
                                    'workspace': params['workspace_name'],
                                    'objects':[{
                                            'type': 'KBaseCollections.FeatureSet',
                                            'data': output_featureSet,
                                            'name': output_name,
                                            'meta': {},
                                            'provenance': provenance
                                        }]
                                })[0]

                    [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
                    objects_created_refs.append(str(new_obj_info[WSID_I])+'/'+str(new_obj_info[OBJID_I]))


            #### Build output report chunks
            ##
            self.log(console,"BUILDING REPORT CHUNK")  # DEBUG
            if len(invalid_msgs) == 0:

                # text report
                #
                report += 'sequences in search db: '+str(seq_total)+"\n"
                report += 'sequences in hit set: '+str(len(hit_order))+"\n"
                report += 'sequences in accepted hit set: '+str(hit_total)+"\n"
                report += "\n"
                #for line in hit_buf:
                #    report += line
                self.log (console, report)


                # build html report chunk
                if many_type_name == 'Genome':
                    feature_id_to_function = GenomeToFASTA_retVal['feature_id_to_function']
                    genome_ref_to_sci_name = GenomeToFASTA_retVal['genome_ref_to_sci_name']
                elif many_type_name == 'GenomeSet':
                    feature_id_to_function = GenomeSetToFASTA_retVal['feature_id_to_function']
                    genome_ref_to_sci_name = GenomeSetToFASTA_retVal['genome_ref_to_sci_name']
                elif many_type_name == 'FeatureSet':
                    feature_id_to_function = FeatureSetToFASTA_retVal['feature_id_to_function']
                    genome_ref_to_sci_name = FeatureSetToFASTA_retVal['genome_ref_to_sci_name']

                head_color = "#eeeeff"
                border_head_color = "#ffccff"
                accept_row_color = 'white'
                #reject_row_color = '#ffeeee'
                reject_row_color = '#eeeeee'
                reject_cell_color = '#ffcccc'
                text_fontsize = "2"
                text_color = '#606060'
                border_body_color = "#cccccc"
                bar_width = 100
                bar_height = 15
                bar_color = "lightblue"
                bar_line_color = "#cccccc"
                bar_fontsize = "1"
                bar_char = "."
                cellpadding = "3"
                cellspacing = "2"
                border = "0"

                html_report_chunk = []

                for line in hit_buf:
                    line = line.strip()
                    if line == '' or line.startswith('#'):
                        continue

                    [hit_id, hit_accession, query_name, query_accession, e_value, bit_score, bias, e_value_best_dom, bit_score_best_dom, bias_best_dom, expected_dom_n, regions, regions_multidom, overlaps, envelopes, dom_n, doms_within_rep_thresh, doms_within_inc_thresh, hit_desc] = re.split('\s+',line)[0:19]

    #                [query_id, hit_id, identity, aln_len, mismatches, gap_openings, q_beg, q_end, h_beg, h_end, e_value, bit_score] = line.split("\t")[0:12]
    #                identity = str(round(float(identity), 1))
    #                if identity == '100.0':  identity = '100'

                    # get coords with respect to hit sequence
                    h_len = hit_seq_len[hit_id]
                    h_beg = hit_beg[hit_id]
                    h_end = hit_end[hit_id]
                    aln_len = h_end-h_beg+1
                    aln_len_perc = round (100.0*float(aln_len)/float(h_len), 1)


                    #if many_type_name == 'SingleEndLibrary':
                    #    pass
                    #elif many_type_name == 'SequenceSet':
                    if many_type_name == 'SequenceSet':
                        pass
                    elif many_type_name == 'Genome' or \
                            many_type_name == 'GenomeSet' or \
                            many_type_name == 'FeatureSet':

                        if many_type_name != 'Genome':
                            [genome_ref, hit_fid] = hit_id.split(genome_id_feature_id_delim)
                        else:
                            genome_ref = input_many_ref
                            hit_fid = hit_id

                        # can't just use hit_fid because may have pipes translated and can't translate back
                        fid_lookup = None
                        for fid in feature_id_to_function[genome_ref].keys():
                            id_untrans = fid
                            id_trans = re.sub ('\|',':',id_untrans)  # BLAST seems to make this translation now when id format has simple 'kb|blah' format

                            #self.log (console, "SCANNING FIDS.  HIT_FID: '"+str(hit_fid)+"' FID: '"+str(fid)+"' TRANS: '"+str(id_trans)+"'")  # DEBUG

                            if id_untrans == hit_fid or id_trans == hit_fid:
                                #self.log (console, "GOT ONE!")  # DEBUG
                                if many_type_name == 'Genome':
                                    accept_id = fid
                                elif many_type_name == 'GenomeSet' or many_type_name == 'FeatureSet':
                                    accept_id = genome_ref+genome_id_feature_id_delim+fid
                                if accept_id in accept_fids:
                                    row_color = accept_row_color
                                else:
                                    row_color = reject_row_color
                                fid_lookup = fid
                                break
                        #self.log (console, "HIT_FID: '"+str(hit_fid)+"' FID_LOOKUP: '"+str(fid_lookup)+"'")  # DEBUG
                        if fid_lookup == None:
                            raise ValueError ("unable to find fid for hit_fid: '"+str(hit_fid))
                        elif fid_lookup not in feature_id_to_function[genome_ref]:
                            raise ValueError ("unable to find function for fid: '"+str(fid_lookup))
                        fid_disp = re.sub (r"^.*\.([^\.]+)\.([^\.]+)$", r"\1.\2", fid_lookup)

                        func_disp = feature_id_to_function[genome_ref][fid_lookup]
                        genome_sci_name = genome_ref_to_sci_name[genome_ref]

                        #if 'overlap_fraction' in params and float(params['overlap_fraction']) > float(high_bitscore_alnlen[hit_seq_id])/float(query_len):

                        html_report_chunk += ['<tr bgcolor="'+row_color+'">']
                        #html_report_chunk += ['<tr bgcolor="'+'white'+'">']  # DEBUG
                        # add overlap bar

                        # coverage graphic (with respect to hit seq)
                        html_report_chunk += ['<td valign=middle align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'">']
                        html_report_chunk += ['<table style="height:'+str(bar_height)+'px; width:'+str(bar_width)+'px" border=0 cellpadding=0 cellspacing=0>']
                        full_len_pos = bar_width
                        aln_beg_pos = int (float(bar_width) * float(int(h_beg)-1)/float(int(h_len)-1))
                        aln_end_pos = int (float(bar_width) * float(int(h_end)-1)/float(int(h_len)-1))
                        cell_pix_height = str(int(round(float(bar_height)/3.0, 0)))

                        cell_color = ['','','']
                        cell_width = []
                        cell_width.append(aln_beg_pos)
                        cell_width.append(aln_end_pos-aln_beg_pos)
                        cell_width.append(bar_width-aln_end_pos)

                        for row_i in range(3):
                            html_report_chunk += ['<tr style="height:'+cell_pix_height+'px">']
                            unalign_color = row_color
                            if row_i == 1:
                                unalign_color = bar_line_color
                            cell_color[0] = unalign_color
                            cell_color[1] = bar_color
                            cell_color[2] = unalign_color

                            for col_i in range(3):
                                cell_pix_width = str(cell_width[col_i])
                                cell_pix_color = cell_color[col_i]
                                html_report_chunk += ['<td style="height:'+cell_pix_height+'px; width:'+cell_pix_width+'px" bgcolor="'+cell_pix_color+'"></td>']
                            html_report_chunk += ['</tr>']
                        html_report_chunk += ['</table>']
                        html_report_chunk += ['</td>']

                        # add other cells
                        # fid
                        html_report_chunk += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(fid_disp)+'</font></td>']
    #                    html_report_chunk += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(hit_accession)+'</font></td>']
                        # func
                        html_report_chunk += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+func_disp+'</font></td>']
                        # sci name
                        html_report_chunk += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+genome_sci_name+'</font></td>']
                        # ident
    #                    if 'ident_thresh' in filtering_fields[hit_id]:
     #                       this_cell_color = reject_cell_color
     #                   else:
     #                       this_cell_color = row_color
     #                   html_report_chunk += ['<td align=center bgcolor="'+this_cell_color+'" style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(identity)+'%</font></td>']

                        # aln len
     #                   if 'overlap_fraction' in filtering_fields[hit_id]:
     #                       this_cell_color = reject_cell_color
     #                   else:
     #                       this_cell_color = row_color
                        html_report_chunk += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(aln_len)+' ('+str(aln_len_perc)+'%)</font></td>']

                        # evalue
                        html_report_chunk += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'><nobr>'+str(e_value)+'</nobr></font></td>']

                        # bit score
                        if 'bitscore' in filtering_fields[hit_id]:
                            this_cell_color = reject_cell_color
                        else:
                            this_cell_color = row_color
                        html_report_chunk += ['<td align=center bgcolor="'+str(this_cell_color)+'" style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'><nobr>'+str(bit_score)+'</nobr></font></td>']
                        # bias
    #                    html_report_chunk += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'><nobr>'+str(bias)+'</nobr><br><nobr>('+str(bias_best_dom)+')</nobr></font></td>']

                        # aln coords only for hit seq
                        html_report_chunk += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'><nobr>'+str(h_beg)+'-'+str(h_end)+'</nobr></font></td>']

                        # close chunk
                        html_report_chunk += ['</tr>']

                # attach chunk
                html_report_chunk_str = "\n".join(html_report_chunk)
                html_report_chunks.append(html_report_chunk_str)
                #self.log(console, "HTML_REPORT_CHUNK: '"+str(html_report_chunk_str)+"'")  # DEBUG


        #### Create and Upload output objects if coalesce_output is true
        ##
        if 'coalesce_output' in params and int(params['coalesce_output']) == 1:
            output_name = params['output_filtered_name']

            if len(invalid_msgs) == 0:
                if not hit_accept_something:
                    self.log(console,"No Object to Upload for all MSAs")  # DEBUG

                else:
                    self.log(console,"Uploading results Object")  # DEBUG

                    if many_type_name == 'SequenceSet':  # input many SequenceSet -> save SequenceSet

                        output_sequenceSet['sequences'] = coalesced_sequenceObjs
                        new_obj_info = ws.save_objects({
                                'workspace': params['workspace_name'],
                                'objects':[{
                                        'type': 'KBaseSequences.SequenceSet',
                                        'data': output_sequenceSet,
                                        'name': output_name,
                                        'meta': {},
                                        'provenance': provenance
                                        }]
                                })[0]

                    else:  # input FeatureSet, Genome, and GenomeSet -> upload FeatureSet output

                        output_featureSet['element_ordering'] = coalesce_featureIds_element_ordering
                        output_featureSet['elements'] = dict()
                        for i,fId in enumerate(output_featureSet['element_ordering']):
                            output_featureSet['elements'][fId] = []
                            output_featureSet['elements'][fId].append(coalesce_featureIds_genome_ordering[i])

                        new_obj_info = ws.save_objects({
                                'workspace': params['workspace_name'],
                                'objects':[{
                                        'type': 'KBaseCollections.FeatureSet',
                                        'data': output_featureSet,
                                        'name': output_name,
                                        'meta': {},
                                        'provenance': provenance
                                        }]
                                })[0]

                    [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
                    objects_created_refs.append(str(new_obj_info[WSID_I])+'/'+str(new_obj_info[OBJID_I]))


        #### Set paths for output HTML
        ##
        html_output_dir = os.path.join(self.output_dir,'html_output')
        if not os.path.exists(html_output_dir):
            os.makedirs(html_output_dir)
        html_search_file = search_tool_name+'_Search.html'
        html_search_path = os.path.join (html_output_dir, html_search_file)
        html_profile_file = search_tool_name+'_Profile.html'
        html_profile_path = os.path.join (html_output_dir, html_profile_file)


        #### Build Search output report (and assemble html chunks)
        ##
        self.log(console,"BUILDING SEARCH REPORT ")  # DEBUG
        if len(invalid_msgs) == 0:

            # build html report
            if many_type_name == 'Genome':
                feature_id_to_function = GenomeToFASTA_retVal['feature_id_to_function']
                genome_ref_to_sci_name = GenomeToFASTA_retVal['genome_ref_to_sci_name']
            elif many_type_name == 'GenomeSet':
                feature_id_to_function = GenomeSetToFASTA_retVal['feature_id_to_function']
                genome_ref_to_sci_name = GenomeSetToFASTA_retVal['genome_ref_to_sci_name']
            elif many_type_name == 'FeatureSet':
                feature_id_to_function = FeatureSetToFASTA_retVal['feature_id_to_function']
                genome_ref_to_sci_name = FeatureSetToFASTA_retVal['genome_ref_to_sci_name']
                
            sp = '&nbsp;'
            head_color = "#eeeeff"
            border_head_color = "#ffccff"
            accept_row_color = 'white'
            #reject_row_color = '#ffeeee'
            reject_row_color = '#eeeeee'
            reject_cell_color = '#ffcccc'
            text_fontsize = "2"
            text_color = '#606060'
            header_tab_fontsize = "3"
            header_tab_color = '#606060'
            border_body_color = "#cccccc"
            bar_width = 100
            bar_height = 15
            bar_color = "lightblue"
            bar_line_color = "#cccccc"
            bar_fontsize = "1"
            bar_char = "."
            cellpadding = "3"
            cellspacing = "2"
            border = "0"
            table_col_width = 8

            html_report_lines = []
            html_report_lines += ['<html>']
            html_report_lines += ['<head>']
            html_report_lines += ['<title>KBase HMMER Custom Model Search Hits</title>']
            html_report_lines += ['</head>']
            html_report_lines += ['<body bgcolor="white">']
            if many_type_name == 'GenomeSet':
                html_report_lines += ['<a href="'+html_profile_file+'" target="profile"><font color="'+header_tab_color+'" size='+header_tab_fontsize+'>TABULAR PROFILE</font></a> | <font color="'+header_tab_color+'" size='+header_tab_fontsize+'><b>SEARCH HITS</b></font>']
                html_report_lines += ['<p>']
            html_report_lines += ['<table cellpadding='+cellpadding+' cellspacing = '+cellspacing+' border='+border+'>']
            html_report_lines += ['<tr bgcolor="'+head_color+'">']
            html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'ALIGNMENT COVERAGE (HIT SEQ)'+'</font></td>']
            html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'GENE ID'+'</font></td>']
            html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'FUNCTION'+'</font></td>']
            html_report_lines += ['<td style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'GENOME'+'</font></td>']
#            html_report_lines += ['<td align=center style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'IDENT'+'%</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'ALN_LEN'+'</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'E-VALUE'+'</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'BIT SCORE'+'</font></td>']
            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'<nobr>H_BEG-H_END</nobr>'+'</font></td>']
#            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'MIS MATCH'+'</font></td>']
#            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'GAP OPEN'+'</font></td>']
            html_report_lines += ['</tr>']


            for msa_i,input_msa_name in enumerate(input_msa_names):
                html_report_lines += ['<tr><td colspan=table_col_width>Hits to <b>'+str(input_msa_name)+'</b></td></tr>']
                if total_hit_cnt[msa_i] == 0:
                    html_report_lines += ['<tr><td colspan=table_col_width><blockquote><i>no hits found</i></td></tr>']
                else:
                    #html_report_lines.extend(html_report_chunks[i])
                    html_report_lines += [ html_report_chunks[msa_i] ]
                html_report_lines += ['<tr><td colspan=table_col_width>'+sp+'</td></tr>']

            html_report_lines += ['</table>']
            html_report_lines += ['</body>']
            html_report_lines += ['</html>']

            # write html to file
            html_path = html_search_path
            html_report_str = "\n".join(html_report_lines)
            with open (html_path, 'w', 0) as html_handle:
                html_handle.write(html_report_str)


        #### Build Profile output report
        ##
        self.log(console,"BUILDING PROFILE REPORT ")  # DEBUG
        if len(invalid_msgs) == 0 and many_type_name == 'GenomeSet':

            # calculate table
            #
            cats = input_msa_names
            table_data = dict()
            INSANE_VALUE = 10000000000000000
            overall_low_val  =  INSANE_VALUE
            overall_high_val = -INSANE_VALUE
            cat_seen = dict()

            # count raw
            for genome_ref in genome_refs:
                if genome_ref not in table_data:
                    table_data[genome_ref] = dict()
                for cat in cats:
                    table_data[genome_ref][cat] = 0
                    
                if genome_ref not in hit_cnt_by_genome_and_model:
                    continue

                for cat in cats:
                    if cat in hit_cnt_by_genome_and_model[genome_ref] and \
                       hit_cnt_by_genome_and_model[genome_ref][cat] != 0:
                        table_data[genome_ref][cat] = hit_cnt_by_genome_and_model[genome_ref][cat]
                        cat_seen[cat] = True

            # determine high and low val
            for genome_ref in genome_refs:
                for cat in cats:
                    val = table_data[genome_ref][cat]
                    if val == 0:  continue
                    #self.log (console, "HIGH VAL SCAN CAT: '"+cat+"' VAL: '"+str(val)+"'")  # DEBUG
                    if val > overall_high_val:
                        overall_high_val = val
                    if val < overall_low_val:
                        overall_low_val = val
            if overall_high_val == -INSANE_VALUE:
                raise ValueError ("unable to find any counts")


            # build html report
            sp = '&nbsp;'
            text_color   = "#606060"
            text_color_2 = "#606060"
            head_color_1 = "#eeeeee"
            head_color_2 = "#eeeeee"
            border_color = "#cccccc"
            border_cat_color = "#ffccff"
            #graph_color = "lightblue"
            #graph_width = 100
            #graph_char = "."
            graph_char = sp
            color_list = ['0','1','2','3','4','5','6','7','8','9','a','b','c','d','e']
            max_color = len(color_list)-1
            cat_disp_trunc_len = 40
            cell_width = '10px'
            if len(genome_refs) > 20:
                graph_gen_fontsize = "1"
#            elif len(genome_refs) > 10:
#                graph_gen_fontsize = "2"
            else:
#                graph_gen_fontsize = "3"
                graph_gen_fontsize = "2"
            if len(cats) > 20:
                graph_cat_fontsize = "1"
#            elif len(cats) > 5:
#                graph_cat_fontsize = "2"
            else:
#                graph_cat_fontsize = "3"
                graph_cat_fontsize = "2"
            if int(graph_cat_fontsize) < int(graph_gen_fontsize):
                cell_fontsize = graph_gen_fontsize = graph_cat_fontsize
            else:
                cell_fontsize = graph_cat_fontsize = graph_gen_fontsize
            graph_padding = "5"
            graph_spacing = "3"
            #border = "1"
            border = "0"
            #row_spacing = "-2"
            num_rows = len(genome_refs)
            show_groups = False
            show_blanks = False
            if 'show_blanks' in params and int(params['show_blanks']) == 1:
                show_blanks = True

            # build html buffer
            html_report_lines = []
            html_report_lines += ['<html>']
            html_report_lines += ['<head>']
            html_report_lines += ['<title>KBase HMMER Custom Model Profile</title>']
            html_report_lines += ['<style>']
            html_report_lines += [".vertical-text {\ndisplay: inline-block;\noverflow: hidden;\nwidth: 0.65em;\n}\n.vertical-text__inner {\ndisplay: inline-block;\nwhite-space: nowrap;\nline-height: 1.1;\ntransform: translate(0,100%) rotate(-90deg);\ntransform-origin: 0 0;\n}\n.vertical-text__inner:after {\ncontent: \"\";\ndisplay: block;\nmargin: 0.0em 0 100%;\n}"]
            html_report_lines += [".vertical-text_title {\ndisplay: inline-block;\noverflow: hidden;\nwidth: 1.0em;\n}\n.vertical-text__inner_title {\ndisplay: inline-block;\nwhite-space: nowrap;\nline-height: 1.0;\ntransform: translate(0,100%) rotate(-90deg);\ntransform-origin: 0 0;\n}\n.vertical-text__inner_title:after {\ncontent: \"\";\ndisplay: block;\nmargin: 0.0em 0 100%;\n}"]
            html_report_lines += ['</style>']
            html_report_lines += ['</head>']
            html_report_lines += ['<body bgcolor="white">']
            html_report_lines += ['<font color="'+header_tab_color+'" size='+header_tab_fontsize+'><b>TABULAR PROFILE</b></font> | <a href="'+html_search_file+'" target="search"><font color="'+header_tab_color+'" size='+header_tab_fontsize+'>SEARCH HITS</font></a>']
            html_report_lines += ['<p>']

            # genomes as rows
            if 'vertical' in params and int(params['vertical']) == 1:
                # table header
                html_report_lines += ['<table cellpadding='+graph_padding+' cellspacing='+graph_spacing+' border='+border+'>']
                corner_rowspan = "1"
                label = ''
                html_report_lines += ['<tr>']
                html_report_lines += ['<td valign=bottom align=right rowspan='+corner_rowspan+'><div class="vertical-text_title"><div class="vertical-text__inner_title"><font color="'+text_color+'">'+label+'</font></div></div></td>']


                # column headers
                for cat in cats:
                    if not cat_seen[cat] and not show_blanks:
                        continue
                    cat_disp = cat
                    cell_title = cat_disp
                    if len(cat_disp) > cat_disp_trunc_len:
                        cat_disp = cat_disp[0:cat_disp_trunc_len]+'*'
                    html_report_lines += ['<td style="border-right:solid 2px '+border_cat_color+'; border-bottom:solid 2px '+border_cat_color+'" bgcolor="'+head_color_2+'"title="'+cell_title+'" valign=bottom align=center>']
                    html_report_lines += ['<div class="vertical-text"><div class="vertical-text__inner">']
                    html_report_lines += ['<font color="'+text_color_2+'" size='+graph_cat_fontsize+'><b>']
                    #for c_i,c in enumerate(cat_disp):
                    #    if c_i < len(cat_disp)-1:
                    #        html_report_lines += [c+'<br>']
                    #    else:
                    #        html_report_lines += [c]
                    html_report_lines += [cat_disp]
                    html_report_lines += ['</b></font>']
                    html_report_lines += ['</div></div>']
                    html_report_lines += ['</td>']
                html_report_lines += ['</tr>']

                # rest of rows
                for genome_ref in genome_refs:
                    genome_sci_name = genome_ref_to_sci_name[genome_ref]
                    html_report_lines += ['<tr>']
                    html_report_lines += ['<td align=right><font color="'+text_color+'" size='+graph_gen_fontsize+'><b><nobr>'+genome_sci_name+'</nobr></b></font></td>']
                    for cat in cats:
                        if not cat_seen[cat] and not show_blanks:
                            continue
                        val = table_data[genome_ref][cat]
                        if val == 0:
                            cell_color = 'white'
                        else:                    
                            cell_color_i = max_color - int(round(max_color * (val-overall_low_val) / float(overall_high_val-overall_low_val)))
                            c = color_list[cell_color_i]
                            cell_color = '#'+c+c+c+c+'FF'

                        cell_val = str(table_data[genome_ref][cat])  # the key line

                        if 'heatmap' in params and params['heatmap'] == '1':
                            if table_data[genome_ref][cat] == 0:
                                this_text_color = text_color
                                #this_graph_char = "0"
                                this_graph_char = sp
                            else:
                                this_text_color = cell_color
                                this_graph_char = graph_char
                                html_report_lines += ['<td align=center valign=middle title="'+cell_val+'" style="width:'+cell_width+'" bgcolor="'+cell_color+'"><font color="'+this_text_color+'" size='+cell_fontsize+'>'+this_graph_char+'</font></td>']
                        else:
                            html_report_lines += ['<td align=center valign=middle style="'+cell_width+'; border-right:solid 2px '+border_color+'; border-bottom:solid 2px '+border_color+'"><font color="'+text_color+'" size='+cell_fontsize+'>'+cell_val+'</font></td>']

                    html_report_lines += ['</tr>']
                html_report_lines += ['</table>']

            # genomes as columns
            else:
                raise ValueError ("Do not yet support Genomes as columns")


            # key table
            html_report_lines += ['<p>']
            html_report_lines += ['<table cellpadding=3 cellspacing=2 border='+border+'>']
            html_report_lines += ['<tr><td valign=middle align=left colspan=3 style="border-bottom:solid 4px '+border_color+'"><font color="'+text_color+'"><b>KEY</b></font></td></tr>']

            for cat_i,cat in enumerate(cats):
                cell_color = 'white'
                if not cat_seen[cat] and not show_blanks:
                    cell_color = "#eeeeee"
                desc = input_msa_descs[cat_i]
                cat_disp = cat
                if len(cat_disp) > cat_disp_trunc_len+1:
                    cat_disp = cat_disp[0:cat_disp_trunc_len]+'*'
                html_report_lines += ['<tr>']
                html_report_lines += ['<td valign=middle align=left bgcolor="'+cell_color+'" style="border-right:solid 4px '+border_color+'><font color="'+text_color+'" size='+graph_cat_fontsize+'>'+cat_disp+'</font></td>']
                html_report_lines += ['<td valign=middle align=left bgcolor="'+cell_color+'"><font color="'+text_color+'" size='+graph_cat_fontsize+'>'+desc+'</font></td>']
                html_report_lines += ['</tr>']

            html_report_lines += ['</table>']

            # close
            html_report_lines += ['</body>']
            html_report_lines += ['</html>']

            # write html to file and upload
            html_path = html_profile_path
            html_report_str = "\n".join(html_report_lines)
            with open (html_path, 'w', 0) as html_handle:
                html_handle.write(html_report_str)


        #### Upload HTML reports
        ##
        self.log(console,"UPLOADING HTML REPORT(s)")  # DEBUG
        if len(invalid_msgs) == 0:

            # Upload HTML Report dir
            #
            dfu = DFUClient(self.callbackURL)
            # upload output html
            try:
                #HTML_upload_ret = dfu.file_to_shock({'file_path': html_path,
                HTML_upload_ret = dfu.file_to_shock({'file_path': html_output_dir,
                                                     'make_handle': 0,
                                                     'pack': 'zip'})
            except:
                raise ValueError ('Logging exception loading HTML file to shock')


        #### Upload output files
        ##
        self.log(console,"UPLOADING OUTPUT FILES")  # DEBUG
        if len(invalid_msgs) == 0:

            output_hit_TAB_dir = os.path.join(self.output_dir, 'HMMER_output_TAB')
            output_hit_MSA_dir = os.path.join(self.output_dir, 'HMMER_output_MSA')
            
            for msa_i,input_msa_name in enumerate(input_msa_names):
                if total_hit_cnt[msa_i] == 0:
                    self.log(console, 'SKIPPING UPLOAD OF EMPTY HMMER OUTPUT FOR MSA '+input_msa_name)
                    continue
                new_hit_TAB_file_path = os.path.join(output_hit_TAB_dir, input_msa_name+'.hitout.txt');
                new_hit_MSA_file_path = os.path.join(output_hit_MSA_dir, input_msa_name+'.msaout.txt');

                shutil.copy(output_hit_TAB_file_paths[msa_i], new_hit_TAB_file_path)
                shutil.copy(output_hit_MSA_file_paths[msa_i], new_hit_MSA_file_path)

            # Upload output dirs
            TAB_upload_ret = None
            MSA_upload_ret = None
            self.log(console, 'UPLOADING OF HMMER OUTPUT FOR MSA '+input_msa_name)
            try:
                TAB_upload_ret = dfu.file_to_shock({'file_path': output_hit_TAB_dir,
                                                    'make_handle': 0,
                                                    'pack': 'zip'})
            except:
                raise ValueError ('Logging exception loading TAB output to shock')
            try:
                MSA_upload_ret = dfu.file_to_shock({'file_path': output_hit_MSA_dir,
                                                    'make_handle': 0,
                                                    'pack': 'zip'})
            except:
                raise ValueError ('Logging exception loading MSA output to shock')


        #### Create report object
        ##
        self.log(console,"CREATING REPORT OBJECT")  # DEBUG
        if len(invalid_msgs) == 0:

            reportName = 'hmmer_report_'+str(uuid.uuid4())
            reportObj = {'objects_created': [],
                         #'text_message': '',  # or is it 'message'?
                         'message': '',  # or is it 'text_message'?
                         'direct_html': None,
                         'direct_html_link_index': None,
                         'file_links': [],
                         'html_links': [],
                         'workspace_name': params['workspace_name'],
                         'report_object_name': reportName
                         }
            #html_buf_lim = 16000  # really 16KB, but whatever
            #if len(html_report_str) <= html_buf_lim:
            #    reportObj['direct_html'] = html_report_str
            #else:

# HERE
            reportObj['direct_html_link_index'] = 0
            reportObj['html_links'] = [{'shock_id': HTML_upload_ret['shock_id'],
                                        'name': html_profile_file,
                                        'label': search_tool_name+' HTML Report'}
                                       ]

            if TAB_upload_ret != None:
                reportObj['file_links'] += [{'shock_id': TAB_upload_ret['shock_id'],
                                             'name': search_tool_name+'_Search.TAB',
                                             'label': search_tool_name+'-'+' hits TABLE'}]
            if MSA_upload_ret != None:
                reportObj['file_links'] += [{'shock_id': MSA_upload_ret['shock_id'],
                                           
                                             'name': search_tool_name+'_Search.MSA',
                                             'label': search_tool_name+' hits MSA'}
                                           ]
            if 'coalesce_output' in params and int(params['coalesce_output']) == 1:
                for object_created_ref in objects_created_refs:
                    reportObj['objects_created'].append({'ref':object_created_ref, 'description':'Coalesced'+' '+search_tool_name+' hits'})
            else:
                for i,input_msa_name in enumerate(input_msa_names):
                    if total_hit_cnt[i] == 0:
                        continue
                    reportObj['objects_created'].append({'ref':objects_created_refs[i], 'description':input_msa_name+' '+search_tool_name+' hits'})
                

            # save report object
            #
            SERVICE_VER = 'release'
            reportClient = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
            #report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})
            report_info = reportClient.create_extended_report(reportObj)


        #### data validation error
        ##
        if len(invalid_msgs) > 0:
            report += "FAILURE\n\n"+"\n".join(invalid_msgs)+"\n"

            reportObj = {
                'objects_created':[],
                'text_message':report
                }

            reportName = 'hmmer_report_'+str(uuid.uuid4())
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
                            'provenance':provenance
                            }
                        ]
                    })[0]
            report_info = dict()
            report_info['name'] = report_obj_info[1]
            report_info['ref'] = str(report_obj_info[6])+'/'+str(report_obj_info[0])+'/'+str(report_obj_info[4])


        #### Return Report
        ##
        self.log(console,"BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = { 'report_name': report_info['name'],
                      'report_ref': report_info['ref']
                      }
        self.log(console,search_tool_name+"_Search DONE")
        #END HMMER_Local_MSA_Group_Search

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method HMMER_Local_MSA_Group_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
