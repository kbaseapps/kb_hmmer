# -*- coding: utf-8 -*-
import json
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

# Tree utils
import ete3

# SDK Utils
from installed_clients.WorkspaceClient import Workspace
from installed_clients.KBaseDataObjectToFileUtilsClient import KBaseDataObjectToFileUtils
from installed_clients.DataFileUtilClient import DataFileUtil as DFUClient
from installed_clients.KBaseReportClient import KBaseReport

# resource Utils
from installed_clients.BBToolsClient import BBTools


# silence whining
import requests
requests.packages.urllib3.disable_warnings()


class HmmerUtil:

    workspaceURL = None
    shockURL = None
    handleURL = None
    callbackURL = None
    scratch = None

    HMMER_BIN = os.path.join(os.sep, 'kb', 'module', 'hmmer', 'bin')
    HMMER_BUILD = os.path.join(HMMER_BIN, 'hmmbuild')  # construct profile HMM(s) from MSA(s)
    HMMER_MAKE_DB = os.path.join(HMMER_BIN, 'makehmmerdb')  # build a HMMER binary db from a seq file
    HMMER_SEARCH = os.path.join(HMMER_BIN, 'hmmsearch')  # search profile(s) against a sequence db

    HMMER_PHMMER = os.path.join(HMMER_BIN, 'phmmer')  # search protein sequence(s) against a protein sequence db
    HMMER_NHMMER = os.path.join(HMMER_BIN, 'nhmmer')  # search nuc sequence(s) against a nuc sequence db
    HMMER_JACKHAMMER = os.path.join(HMMER_BIN, 'jackhmmer')  # iteratively search sequence(s) against a protein db

    #HMMER_ALIGN      = '/kb/module/hmmer/binaries/hmmalign'  # align sequences to a profile HMM
    #HMMER_PRESS      = '/kb/module/hmmer/binaries/hmmpress'  # prepare HMM db for hmmscan
    #HMMER_SCAN       = '/kb/module/hmmer/binaries/hmmscan'  # scan prot sequence(s) against protein profile db
    #HMMER_NSCAN       = '/kb/module/hmmer/binaries/nhmmscan'  # scan nuc sequence(s) against nuc profile db


    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config, ctx):
        self.config = config
        self.ctx = ctx
        
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.serviceWizardURL = config['service-wizard-url']

#        self.callbackURL = os.environ['SDK_CALLBACK_URL'] if os.environ['SDK_CALLBACK_URL'] != None else 'https://kbase.us/services/njs_wrapper'
        self.callbackURL = config['SDK_CALLBACK_URL']
        if self.callbackURL == None:
            raise ValueError("SDK_CALLBACK_URL not set in environment")

        try:
            self.wsClient = Workspace(self.workspaceURL, token=ctx['token'])
        except:
            raise ValueError ("unable to connect to Workspace service")
            
        self.scratch = os.path.abspath(config['scratch'])
        if self.scratch == None:
            self.scratch = os.path.join('/kb', 'module', 'local_scratch')
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        # set i/o dirs
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds() * 1000)
        self.input_dir = os.path.join(self.scratch, 'input.' + str(timestamp))
        self.output_dir = os.path.join(self.scratch, 'output.' + str(timestamp))
        if not os.path.exists(self.input_dir):
            os.makedirs(self.input_dir)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)


    # timestamp
    def now_ISO(self):
        now_timestamp = datetime.now()
        now_secs_from_epoch = (now_timestamp - datetime(1970,1,1)).total_seconds()
        now_timestamp_in_iso = datetime.fromtimestamp(int(now_secs_from_epoch)).strftime('%Y-%m-%d_%T')
        return now_timestamp_in_iso

    # message logging
    def log(self, target, message):
        message = '['+self.now_ISO()+'] '+message
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    # set provenance
    def _set_provenance (self, ctx, service, method, input_obj_refs):
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].extend(input_obj_refs)
        provenance[0]['service'] = service
        provenance[0]['method'] = method
        return provenance
    
    # error report
    def _create_error_report (self, workspace_name, message, provenance):
        reportObj = {
            'objects_created': [],
            'text_message': message
        }

        reportName = 'hmmer_report_' + str(uuid.uuid4())
        report_obj_info = self.wsClient.save_objects({
            #                'id':info[6],
            'workspace': workspace_name,
            'objects': [
                {
                    'type': 'KBaseReport.Report',
                    'data': reportObj,
                    'name': reportName,
                    'meta': {},
                    'hidden': 1,
                    'provenance': provenance
                }
            ]
        })[0]
        returnVal = {'report_name': report_obj_info[1],
                     'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4])
                     }
        return returnVal
        
    # split genome from feature id if id from set
    def _parse_genome_and_feature_id_from_hit_id(self,
                                                 hit_id,
                                                 target_type,
                                                 target_ref,
                                                 genome_id_feature_id_delim):
        genome_ref = None
        feature_id = None
        if target_type == 'Genome' or target_type == 'AnnotatedMetagenomeAssembly':
            genome_ref = target_ref
            feature_id = hit_id
        else:
            [genome_ref, feature_id] = hit_id.split(genome_id_feature_id_delim)
        return [genome_ref, feature_id]

    # read model group version
    def _get_version(self, this_model_group, this_model_group_data_dir):
        this_version = None
        this_version_file = os.path.join(this_model_group_data_dir,this_model_group+'.version')
        if not os.path.exists(this_version_file):
            raise ValueError ("ABORT: missing version file "+this_version_file)
        with open (this_version_file, 'r') as version_handle:
            for version_line in version_handle.readlines():
                version_line = version_line.strip()
                if len(version_line) == 0:
                    continue
                if version_line.startswith('#'):
                    continue
                this_version = version_line
                break
        if this_version == None:
            raise ValueError ("ABORT: No version found in file "+this_version_file)
        return this_version
    

    # read model group disp name
    def _get_disp_name(self, this_model_group, this_model_group_data_dir):
        this_disp_name = None
        this_disp_name_file = os.path.join(this_model_group_data_dir,this_model_group+'.display_name')
        if not os.path.exists(this_disp_name_file):
            raise ValueError ("ABORT: missing display name file "+this_disp_name_file)
        with open (this_disp_name_file, 'r') as disp_name_handle:
            for disp_name_line in disp_name_handle.readlines():
                disp_name_line = disp_name_line.strip()
                if len(disp_name_line) == 0:
                    continue
                if disp_name_line.startswith('#'):
                    continue
                this_disp_name = disp_name_line
                break
        if this_disp_name == None:
            raise ValueError ("ABORT: No display name found in file "+this_disp_name_file)
        return this_disp_name
    

    # read model group's fam groups
    def _get_fam_groups(self, this_model_group, this_hmms_dir):
        this_fam_groups = []
        this_fam_groups_disp = dict()
        this_fam_groups_file = os.path.join(this_hmms_dir,this_model_group+'-categories.txt')
        if not os.path.exists(this_fam_groups_file):
            raise ValueError ("ABORT: missing categories file "+this_fam_groups_file)
        with open (this_fam_groups_file, 'r') as fam_groups_handle:
            for fam_group_line in fam_groups_handle.readlines():
                fam_group_line = fam_group_line.strip()
                if len(fam_group_line) == 0:
                    continue
                if fam_group_line.startswith('#'):
                    continue
                this_fam_group_info = fam_group_line.split("\t")
                this_fam_group = this_fam_group_info[0]
                this_fam_groups.append(this_fam_group)
                if len(this_fam_group_info) > 1:
                    this_fam_groups_disp[this_fam_group] = this_fam_group_info[1]
                else:
                    this_fam_groups_disp[this_fam_group] = this_fam_group_info[0]
                    
        if len(this_fam_groups) == 0:
            raise ValueError ("ABORT: No fam groups found in "+this_fam_groups_file)
        return (this_fam_groups, this_fam_groups_disp)
    

    # set model group config
    def _configure_by_model_group (self, this_model_group):
        if this_model_group == None or this_model_group == '':
            raise ValueError ('ABORT: missing model_group')
        this_cfg = dict()
        this_model_group_data_dir = os.path.join(os.sep, 'kb', 'module', 'data', this_model_group)
        this_version = self._get_version(this_model_group, this_model_group_data_dir)
        this_model_group_disp_name = self._get_disp_name(this_model_group, this_model_group_data_dir)
        this_hmms_dir = os.path.join(this_model_group_data_dir, this_model_group+'-'+this_version)
        this_hmms_path = os.path.join(this_hmms_dir, this_model_group+'-fam-HMMs-'+this_version+'.txt')
        this_cfg = { 'search_tool_name': 'HMMER_'+this_model_group,
                     'group_name': this_model_group_disp_name,
                     'version': this_version,
                     'HMMS_DIR': this_hmms_dir,
                     'HMMS_PATH': this_hmms_path
                   }
        (this_cfg['fam_groups'], this_cfg['fam_groups_disp']) = self._get_fam_groups(this_model_group, this_hmms_dir)
        return this_cfg
    

    # validate MSA type
    def _check_MSA_sequence_type_correct(self, MSA_in, row_order, seq_type):
        PROT_MSA_pattern = re.compile("^[\.\-_acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
        DNA_MSA_pattern = re.compile("^[\.\-_ACGTUXNRYSWKMBDHVacgtuxnryswkmbdhv \t\n]+$")
        this_appropriate_sequence_found_in_MSA_input = True
        msa_invalid_msgs = []

        # Check for PROTEIN sequence type
        #
        if seq_type.startswith('P') or seq_type.startswith('p'):
            if 'sequence_type' in MSA_in and (MSA_in['sequence_type'] == 'dna' or MSA_in['sequence_type'] == 'DNA'):
                this_appropriate_sequence_found_in_MSA_input = False
            else:
                for row_id in row_order:
                    #self.log(console, row_id+": '"+MSA_in['alignment'][row_id]+"'")    # DEBUG
                    if DNA_MSA_pattern.match(MSA_in['alignment'][row_id]):
                        self.log(msa_invalid_msgs,
                                 "Finding nucleotide instead of protein sequences in MSA. " +
                                 "BAD record for MSA row_id: " + row_id + "\n" + MSA_in['alignment'][row_id] + "\n")
                        this_appropriate_sequence_found_in_MSA_input = False
                        break
                    elif not PROT_MSA_pattern.match(MSA_in['alignment'][row_id]):
                        self.log(msa_invalid_msgs,
                                 "Not finding protein sequence in MSA. " +
                                 "BAD record for MSA row_id: " + row_id + "\n" + MSA_in['alignment'][row_id] + "\n")
                        this_appropriate_sequence_found_in_MSA_input = False
                        break

        # Check for NUCLEOTIDE sequence type
        #
        elif seq_type.startswith('N') or seq_type.startswith('n'):
            if 'sequence_type' in MSA_in and (MSA_in['sequence_type'] != 'dna' and MSA_in['sequence_type'] != 'DNA'):
                this_appropriate_sequence_found_in_MSA_input = False
            else:
                for row_id in row_order:
                    #self.log(console, row_id+": '"+MSA_in['alignment'][row_id]+"'")    # DEBUG
                    if not DNA_MSA_pattern.match(MSA_in['alignment'][row_id]):
                        self.log(msa_invalid_msgs,
                                 "Not Finding nucleotide in MSA. " +
                                 "BAD record for MSA row_id: " + row_id + "\n" + MSA_in['alignment'][row_id] + "\n")
                        this_appropriate_sequence_found_in_MSA_input = False
                        break
                    elif PROT_MSA_pattern.match(MSA_in['alignment'][row_id]):
                        self.log(msa_invalid_msgs,
                                 "Finding protein sequence instead of nucleotide sequences in MSA. " +
                                 "BAD record for MSA row_id: " + row_id + "\n" + MSA_in['alignment'][row_id] + "\n")
                        this_appropriate_sequence_found_in_MSA_input = False
                        break

        else:
            raise ValueError("Incorrectly formatted call of _check_MSA_sequence_type_correct() method")

        # return sequence type check logical
        #
        return (this_appropriate_sequence_found_in_MSA_input, msa_invalid_msgs)


    def _get_genome_disp_name (self, 
                               params, 
                               genome_ref, 
                               many_type_name, 
                               ama_ref_to_obj_name, 
                               genome_ref_to_obj_name,
                               genome_ref_to_sci_name):
        if many_type_name == 'AnnotatedMetagenomeAssembly':
            genome_disp_name = ama_ref_to_obj_name[genome_ref]
        else:
            genome_obj_name = genome_ref_to_obj_name[genome_ref]
            genome_sci_name = genome_ref_to_sci_name[genome_ref]
            [ws_id, obj_id, genome_obj_version] = genome_ref.split('/')
            genome_disp_name = ''
            if 'obj_name' in params['genome_disp_name_config']:
                genome_disp_name += genome_obj_name
            if 'ver' in params['genome_disp_name_config']:
                genome_disp_name += '.v'+str(genome_obj_version)
            if 'sci_name' in params['genome_disp_name_config']:
                genome_disp_name += ': '+genome_sci_name

        return genome_disp_name


    def run_HMMER_Model_Group_Search(self, params):
        """
        Method for HMMER search of a meaningful group of Hidden Markov Models
        """
        ctx = self.ctx
        console = []
        invalid_msgs = []
        msa_invalid_msgs = []

        model_group_config = self._configure_by_model_group(params['model_group'])
        
        search_tool_name = model_group_config['search_tool_name']
        self.log(console, 'Running ' + search_tool_name + '_Search with params=')
        self.log(console, "\n" + pformat(params))
        report = ''
#        report = 'Running '+search_tool_name+'_Search with params='
#        report += "\n"+pformat(params)
        #appropriate_sequence_found_in_one_input = False
        #appropriate_sequence_found_in_MSA_input = False
        genome_id_feature_id_delim = '.f:'

        #### do some basic checks
        #
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
#        if 'input_one_ref' not in params:
#            raise ValueError('input_one_ref parameter is required')
#        if 'input_msa_refs' not in params or len(params['input_msa_refs']) == 0:
#            raise ValueError('input_msa_refs parameter is required if selecting local MSAs')
        if 'input_many_refs' not in params:
            raise ValueError('input_many_refs parameter is required')
        if 'genome_disp_name_config' not in params:
            raise ValueError('genome_disp_name_config parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')
        #if 'coalesce_output' not in params:
        #    raise ValueError('coalesce_output parameter is required')

        # set local names and ids
#        input_one_ref = params['input_one_ref']
        #input_msa_ref = params['input_msa_ref']

        input_many_refs = params['input_many_refs']
        ws_id = input_many_refs[0].split('/')[0]

        # set provenance
        self.log(console, "SETTING PROVENANCE")  # DEBUG
        service = 'kb_hmmer'
        method = search_tool_name+'_Search'
        input_obj_refs = params['input_many_refs']
        provenance = self._set_provenance (ctx, service, method, input_obj_refs)
        
        
        # validate that at least one fam is selected and store which fam_ids are explicitly config
        explicitly_requested_models = dict()
        fam_groups = model_group_config['fam_groups']
        some_fam_found = False
        for fam_group in fam_groups:
            fam_field = 'input_'+params['model_group']+'_'+fam_group+'_ids'
            if not params.get(fam_field):
                params[fam_field] = ['ALL']
            for fam_id in params[fam_field]:
                if fam_id.upper() == 'NONE':
                    continue
                if fam_id.upper() != 'ALL':
                    self.log(console, 'RECORDING EXPLICITLY REQUESTED MODEL '+fam_id)  # DEBUG
                    explicitly_requested_models[fam_id] = True                        
                some_fam_found = True
        if not some_fam_found:
            message = 'You must request at least one HMM'
            self.log(invalid_msgs, message)
            return self._create_error_report (params['workspace_name'], message, provenance)


        #### Get the input_many objects
        ##
        all_genome_refs = []
        genome_refs = dict()
        feature_ids = dict()
        genome_CDS_count_by_ref = dict()
        appropriate_sequence_found_in_many_inputs = dict()
        input_many_names = dict()
        many_type_names = dict()
        input_many_objs = dict()
        many_forward_reads_file_paths = dict()
        many_forward_reads_file_handles = dict()
        obj2file_retVal = dict()
        
        for input_many_ref in input_many_refs:

            appropriate_sequence_found_in_many_inputs[input_many_ref] = False

            try:
                #objects = ws.get_objects([{'ref': input_many_ref}])
                objects = self.wsClient.get_objects2({'objects': [{'ref': input_many_ref}]})['data']
                input_many_data = objects[0]['data']
                info = objects[0]['info']
                input_many_name = str(info[1])
                many_type_name = info[2].split('.')[1].split('-')[0]
                input_many_names[input_many_ref] = input_many_name
                many_type_names[input_many_ref] = many_type_name
                
            except Exception as e:
                raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()

            # Handle overloading (input_many can be SequenceSet, FeatureSet, Genome, GenomeSet, Tree, or AMA)
            #
            if many_type_name == 'SequenceSet':
                try:
                    input_many_sequenceSet = input_many_data
                except Exception as e:
                    print(traceback.format_exc())
                    raise ValueError('Unable to get SequenceSet: ' + str(e))

                header_id = input_many_sequenceSet['sequences'][0]['sequence_id']
                many_forward_reads_file_paths[input_many_ref] = os.path.join(self.output_dir, header_id + '.fasta')
                many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
                self.log(console, 'writing reads file: ' + str(many_forward_reads_file_path))

                for seq_obj in input_many_sequenceSet['sequences']:
                    header_id = seq_obj['sequence_id']
                    sequence_str = seq_obj['sequence']

                    PROT_pattern = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
                    DNA_pattern = re.compile("^[acgtuACGTUnryNRY ]+$")
                    if DNA_pattern.match(sequence_str):
                        self.log(invalid_msgs,
                                 "Require protein sequences for target. " +
                                 "BAD nucleotide record for sequence_id: " + header_id + "\n" + sequence_str + "\n")
                        continue
                    elif not PROT_pattern.match(sequence_str):
                        self.log(invalid_msgs, "BAD record for sequence_id: " + header_id + "\n" + sequence_str + "\n")
                        continue
                    appropriate_sequence_found_in_many_inputs[input_many_ref] = True
                    many_forward_reads_file_handle.write('>' + header_id + "\n")
                    many_forward_reads_file_handle.write(sequence_str + "\n")
                many_forward_reads_file_handle.close()
                self.log(console, 'done')


            # we're going to profile genome refs for all target types, even if only one object.
            #  note: we are calling it 'genome_refs' even if the object is an AMA
            #
            genome_refs[input_many_ref] = []

            
            # FeatureSet
            #
            if many_type_name == 'FeatureSet':
                # retrieve sequences for features
                input_many_featureSet = input_many_data
                input_many_objs[input_many_ref] = input_many_featureSet
                many_forward_reads_file_dir = self.output_dir
                many_forward_reads_file = input_many_names[input_many_ref] + ".fasta"

                # DEBUG
                #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
                FeatureSetToFASTA_params = {
                    'featureSet_ref': input_many_ref,
                    'file': many_forward_reads_file,
                    'dir': many_forward_reads_file_dir,
                    'console': console,
                    'invalid_msgs': invalid_msgs,
                    'residue_type': 'protein',
                    'feature_type': 'CDS',
                    'record_id_pattern': '%%genome_ref%%' + genome_id_feature_id_delim + '%%feature_id%%',
                    'record_desc_pattern': '[%%genome_ref%%]',
                    'case': 'upper',
                    'linewrap': 50,
                    'merge_fasta_files': 'TRUE'
                }

                #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
                SERVICE_VER = 'release'
                #SERVICE_VER = 'dev'
                DOTFU = KBaseDataObjectToFileUtils(url=self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
                FeatureSetToFASTA_retVal = DOTFU.FeatureSetToFASTA(FeatureSetToFASTA_params)
                obj2file_retVal[input_many_ref] = FeatureSetToFASTA_retVal
                many_forward_reads_file_paths[input_many_ref] = FeatureSetToFASTA_retVal['fasta_file_path']
                feature_ids[input_many_ref] = dict()
                feature_ids[input_many_ref]['by_genome_ref'] = FeatureSetToFASTA_retVal['feature_ids_by_genome_ref']
                if len(feature_ids[input_many_ref]['by_genome_ref'].keys()) > 0:
                    appropriate_sequence_found_in_many_inputs[input_many_ref] = True
                    these_genome_refs = sorted(feature_ids[input_many_ref]['by_genome_ref'].keys())
                    genome_refs[input_many_ref] = these_genome_refs
                    all_genome_refs.extend(these_genome_refs)
                    
                    for this_genome_ref in these_genome_refs:
                        genome_CDS_count_by_ref[this_genome_ref] = len(feature_ids[input_many_ref]['by_genome_ref'][this_genome_ref])

                # DEBUG
                #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
                #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")

            # Genome
            #
            elif many_type_name == 'Genome':
                many_forward_reads_file_dir = self.output_dir
                many_forward_reads_file = input_many_name + ".fasta"

                # DEBUG
                #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
                GenomeToFASTA_params = {
                    'genome_ref': input_many_ref,
                    'file': many_forward_reads_file,
                    'dir': many_forward_reads_file_dir,
                    'console': console,
                    'invalid_msgs': invalid_msgs,
                    'residue_type': 'protein',
                    'feature_type': 'CDS',
                    'record_id_pattern': '%%feature_id%%',
                    'record_desc_pattern': '[%%genome_id%%]',
                    'case': 'upper',
                    'linewrap': 50
                }

                #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
                #SERVICE_VER = 'release'
                SERVICE_VER = 'dev'
                DOTFU = KBaseDataObjectToFileUtils(url=self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
                GenomeToFASTA_retVal = DOTFU.GenomeToFASTA(GenomeToFASTA_params)
                obj2file_retVal[input_many_ref] = GenomeToFASTA_retVal
                many_forward_reads_file_paths[input_many_ref] = GenomeToFASTA_retVal['fasta_file_path']
                feature_ids[input_many_ref] = dict()
                feature_ids[input_many_ref]['basic'] = GenomeToFASTA_retVal['feature_ids']
                if len(feature_ids[input_many_ref]['basic']) > 0:
                    appropriate_sequence_found_in_many_inputs[input_many_ref] = True
                    genome_refs[input_many_ref] = [input_many_ref]
                    all_genome_refs.append(input_many_ref)
                    
                    genome_CDS_count_by_ref[input_many_ref] = len(feature_ids[input_many_ref]['basic'])

                # DEBUG
                #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
                #self.log(console, "Genome2Fasta() took "+str(end_time-beg_time)+" secs")

                
            # GenomeSet
            #
            elif many_type_name == 'GenomeSet':
                input_many_genomeSet = input_many_data
                input_many_objs[input_many_ref] = input_many_genomeSet
                many_forward_reads_file_dir = self.output_dir
                many_forward_reads_file = input_many_name + ".fasta"

                for genome_id in input_many_genomeSet['elements']:
                    genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                    if genome_ref not in genome_refs:
                        genome_refs[input_many_ref].append(genome_ref)
                        all_genome_refs.append(genome_ref)
                        
                        # DEBUG
                #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
                GenomeSetToFASTA_params = {
                    'genomeSet_ref': input_many_ref,
                    'file': many_forward_reads_file,
                    'dir': many_forward_reads_file_dir,
                    'console': console,
                    'invalid_msgs': invalid_msgs,
                    'residue_type': 'protein',
                    'feature_type': 'CDS',
                    'record_id_pattern': '%%genome_ref%%' + genome_id_feature_id_delim + '%%feature_id%%',
                    'record_desc_pattern': '[%%genome_ref%%]',
                    'case': 'upper',
                    'linewrap': 50,
                    'merge_fasta_files': 'TRUE'
                }

                #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
                #SERVICE_VER = 'release'
                SERVICE_VER = 'dev'
                DOTFU = KBaseDataObjectToFileUtils(url=self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
                GenomeSetToFASTA_retVal = DOTFU.GenomeSetToFASTA(GenomeSetToFASTA_params)
                obj2file_retVal[input_many_ref] = GenomeSetToFASTA_retVal
                many_forward_reads_file_paths[input_many_ref] = GenomeSetToFASTA_retVal['fasta_file_path_list'][0]
                feature_ids[input_many_ref] = dict()
                feature_ids[input_many_ref]['by_genome_id'] = GenomeSetToFASTA_retVal['feature_ids_by_genome_id']
                if len(feature_ids[input_many_ref]['by_genome_id'].keys()) > 0:
                    appropriate_sequence_found_in_many_inputs[input_many_ref] = True

                    for this_genome_id in feature_ids[input_many_ref]['by_genome_id'].keys():
                        this_genome_ref = input_many_genomeSet['elements'][this_genome_id]['ref']
                        genome_CDS_count_by_ref[this_genome_ref] = len(feature_ids[input_many_ref]['by_genome_id'][this_genome_id])

                # DEBUG
                #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
                #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")

                
            # SpeciesTree
            #
            elif many_type_name == 'Tree':
                input_many_tree = input_many_data
                input_many_objs[input_many_ref] = input_many_tree
                many_forward_reads_file_dir = self.output_dir
                many_forward_reads_file = input_many_name + ".fasta"

                if 'type' not in input_many_tree or \
                   input_many_tree['type'] != 'SpeciesTree':
                    raise ValueError ("Tree type must be SpeciesTree for object "+input_many_names[input_many_ref])
                if 'ws_refs' not in input_many_tree:
                    raise ValueError ("Tree must contain ws_refs data for object "+input_many_names[input_many_ref])

                for genome_id in input_many_tree['ws_refs']:
                    genome_ref = input_many_tree['ws_refs'][genome_id]['g'][0]
                    if genome_ref not in genome_refs:
                        genome_refs[input_many_ref].append(genome_ref)
                        all_genome_refs.append(genome_ref)
                        
                # DEBUG
                #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
                SpeciesTreeToFASTA_params = {
                    'tree_ref': input_many_ref,
                    'file': many_forward_reads_file,
                    'dir': many_forward_reads_file_dir,
                    'console': console,
                    'invalid_msgs': invalid_msgs,
                    'residue_type': 'protein',
                    'feature_type': 'CDS',
                    'record_id_pattern': '%%genome_ref%%' + genome_id_feature_id_delim + '%%feature_id%%',
                    'record_desc_pattern': '[%%genome_ref%%]',
                    'case': 'upper',
                    'linewrap': 50,
                    'merge_fasta_files': 'TRUE'
                }

                #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
                #SERVICE_VER = 'release'
                SERVICE_VER = 'dev'
                DOTFU = KBaseDataObjectToFileUtils(url=self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
                SpeciesTreeToFASTA_retVal = DOTFU.SpeciesTreeToFASTA(SpeciesTreeToFASTA_params)
                obj2file_retVal[input_many_ref] = SpeciesTreeToFASTA_retVal
                many_forward_reads_file_paths[input_many_ref] = SpeciesTreeToFASTA_retVal['fasta_file_path_list'][0]
                feature_ids[input_many_ref] = dict()
                feature_ids[input_many_ref]['by_genome_id'] = SpeciesTreeToFASTA_retVal['feature_ids_by_genome_id']
                if len(feature_ids[input_many_ref]['by_genome_id'].keys()) > 0:
                    appropriate_sequence_found_in_many_inputs[input_many_ref] = True

                    for this_genome_id in feature_ids[input_many_ref]['by_genome_id'].keys():
                        this_genome_ref = input_many_tree['ws_refs'][this_genome_id]['g'][0]
                        genome_CDS_count_by_ref[this_genome_ref] = len(feature_ids[input_many_ref]['by_genome_id'][this_genome_id])

                # DEBUG
                #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
                #self.log(console, "FeatureSetToFasta() took "+str(end_time-beg_time)+" secs")

                
            # AnnotatedMetagenomeAssembly
            #
            elif many_type_name == 'AnnotatedMetagenomeAssembly':
                many_forward_reads_file_dir = self.output_dir
                many_forward_reads_file = input_many_name + ".fasta"

                # DEBUG
                #beg_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
                AnnotatedMetagenomeAssemblyToFASTA_params = {
                    'ama_ref': input_many_ref,
                    'file': many_forward_reads_file,
                    'dir': many_forward_reads_file_dir,
                    'console': console,
                    'invalid_msgs': invalid_msgs,
                    'residue_type': 'protein',
                    'feature_type': 'CDS',
                    'record_id_pattern': '%%feature_id%%',
                    'record_desc_pattern': '[%%genome_id%%]',
                    'case': 'upper',
                    'linewrap': 50
                }

                #self.log(console,"callbackURL='"+self.callbackURL+"'")  # DEBUG
                #SERVICE_VER = 'release'
                SERVICE_VER = 'beta'
                DOTFU = KBaseDataObjectToFileUtils(url=self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
                AnnotatedMetagenomeAssemblyToFASTA_retVal = DOTFU.AnnotatedMetagenomeAssemblyToFASTA (AnnotatedMetagenomeAssemblyToFASTA_params)
                obj2file_retVal[input_many_ref] = AnnotatedMetagenomeAssemblyToFASTA_retVal
                many_forward_reads_file_paths[input_many_ref] = AnnotatedMetagenomeAssemblyToFASTA_retVal['fasta_file_path']
                feature_ids[input_many_ref] = dict()
                feature_ids[input_many_ref]['basic'] = AnnotatedMetagenomeAssemblyToFASTA_retVal['feature_ids']
                if len(feature_ids) > 0:
                    appropriate_sequence_found_in_many_inputs[input_many_ref] = True
                    genome_CDS_count_by_ref[input_many_ref] = len(feature_ids[input_many_ref]['basic'])

                genome_refs[input_many_ref] = [input_many_ref]
                all_genome_refs.append(input_many_ref)
                
                # DEBUG
                #end_time = (datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()
                #self.log(console, "Genome2Fasta() took "+str(end_time-beg_time)+" secs")

            # Missing proper input_many_type
            #
            else:
                raise ValueError('Cannot yet handle input_many type of: ' + many_type_name)

        # Get total number of sequences in input_many search db
        #
        seq_total = 0
        for input_many_ref in input_many_refs:
            if many_type_names[input_many_ref] == 'SequenceSet':
                seq_total = len(input_many_sequenceSet['sequences'])
            elif many_type_names[input_many_ref] == 'FeatureSet':
                seq_total = len(input_many_featureSet['elements'].keys())
            elif many_type_names[input_many_ref] == 'Genome' \
                 or many_type_names[input_many_ref] == 'AnnotatedMetagenomeAssembly':
                seq_total = len(feature_ids[input_many_ref]['basic'])
            elif many_type_names[input_many_ref] == 'GenomeSet':
                for genome_id in feature_ids[input_many_ref]['by_genome_id'].keys():
                    seq_total += len(feature_ids[input_many_ref]['by_genome_id'][genome_id])
            elif many_type_names[input_many_ref] == 'Tree':
                for genome_id in feature_ids[input_many_ref]['by_genome_id'].keys():
                    seq_total += len(feature_ids[input_many_ref]['by_genome_id'][genome_id])

                
        #### Extract the HMMs into buf
        ##
        HMM_bufs = dict()
        this_buf = []
        this_id = None
        HMMS_PATH = model_group_config['HMMS_PATH']
        with open(HMMS_PATH, 'r', 0) as hmm_handle:
            for hmm_line in hmm_handle.readlines():
                hmm_line = hmm_line.rstrip()
                if hmm_line.startswith('//'):
                    if len(this_buf) > 0:
                        if this_id == None:
                            raise ValueError("Failure parsing file: " + HMMS_PATH)
                        this_buf.append(hmm_line)
                        HMM_bufs[this_id] = this_buf
                        this_buf = []
                        this_id = None
                        continue
                if hmm_line.startswith('NAME'):
                    this_id = hmm_line.split()[1].replace('.hmm', '')
                this_buf.append(hmm_line)

                
        #### Get all the HMM cats and ids
        ##
        all_HMM_groups_order = []
        all_HMM_ids = dict()
        all_HMM_ids_order = []
        input_HMM_descs = dict()
        model_len = dict()
        hmm_group_config_path = os.path.join(model_group_config['HMMS_DIR'], params['model_group']+'-categories.txt')
        HMM_fam_config_dir = os.path.join(model_group_config['HMMS_DIR'], params['model_group']+'-fams')
        HMM_fam_input_dir = os.path.join(self.output_dir, 'HMMs')

        # get connections between hmms and groups
        with open(hmm_group_config_path, 'r', 0) as hmm_group_config_handle:
            for hmm_group_config_line in hmm_group_config_handle.readlines():
                hmm_group_config_line = hmm_group_config_line.rstrip()
                hmm_group = hmm_group_config_line.split("\t")[0]
                all_HMM_groups_order.append(hmm_group)
                all_HMM_ids[hmm_group] = []
                #HMM_fam_config_path = os.path.join(HMM_fam_config_dir, 'dbCAN-' + hmm_group + '.txt')
                HMM_fam_config_path = os.path.join(HMM_fam_config_dir, params['model_group']+'-' + hmm_group + '.desc.tsv')
                with open(HMM_fam_config_path, 'r', 0) as hmm_fam_config_handle:
                    for hmm_fam_config_line in hmm_fam_config_handle.readlines():
                        hmm_fam_config_line = hmm_fam_config_line.rstrip()
                        hmm_fam_config = hmm_fam_config_line.split("\t")
                        hmm_fam_id = hmm_fam_config[0]
                        if len(hmm_fam_config) > 1:
                            input_HMM_descs[hmm_fam_id] = re.sub(r'[^\x00-\x7F]+',' ', hmm_fam_config[1])
                        else:
                            input_HMM_descs[hmm_fam_id] = hmm_fam_id
                        all_HMM_ids[hmm_group].append(hmm_fam_id)
                        all_HMM_ids_order.append(hmm_fam_id)
                        #self.log(console,"HMM_FAM_CONFIG: "+hmm_group+" "+hmm_fam_id) # DEBUG
                        
        #### get the specific input HMM ids requested
        ##
        input_HMM_ids = dict()
        for hmm_group in all_HMM_groups_order:
            input_HMM_ids[hmm_group] = []
            input_field = 'input_'+params['model_group']+'_' + hmm_group + '_ids'
            if input_field in params and params[input_field] != None and len(params[input_field]) > 0:
                only_none_found = True
                for HMM_fam in params[input_field]:
                    if HMM_fam.upper() == 'NONE':
                        continue
                    elif HMM_fam.upper() == 'ALL':
                        input_HMM_ids[hmm_group] = all_HMM_ids[hmm_group]
                        only_none_found = False
                        break
                    else:
                        only_none_found = False
                        input_HMM_ids[hmm_group].append(HMM_fam)
                if only_none_found:
                    self.log(console, "No HMMs configured for GROUP: "+hmm_group)  # DEBUG
                    input_HMM_ids[hmm_group] = []
            else:  # default: use NONE
                self.log(console, "SKIPPING GROUP: "+hmm_group)  # DEBUG
                input_HMM_ids[hmm_group] = []
            # DEBUG
            #for HMM_fam in input_HMM_ids[hmm_group]:
            #    self.log(console, "SEARCHING WITH "+hmm_group+" "+HMM_fam)

        # check for failed input file creation
        #
        for input_many_ref in input_many_refs:
            if not appropriate_sequence_found_in_many_inputs[input_many_ref]:
                self.log(invalid_msgs, "no protein sequences found in '" + input_many_names[input_many_ref] + "'")

        # input data failed validation.  Need to return
        #
        if len(invalid_msgs) > 0:
            message = "\n".join(invalid_msgs)
            return self._create_error_report (params['workspace_name'], message, provenance)

        
        #### Iterate through categories and make separate Search HITs for each category
        ##
        hmm_groups_used = []
        for hmm_group in all_HMM_groups_order:
            if hmm_group not in input_HMM_ids or input_HMM_ids[hmm_group] == None or len(input_HMM_ids[hmm_group]) == 0:
                continue
            else:
                hmm_groups_used.append(hmm_group)

        # get mapping to from hmm_id to hmm_groups used (for link to search html file)
        hmm_id_to_hmm_group_name = dict()
        hmm_id_to_hmm_group_index = dict()
        for hmm_group_i, hmm_group in enumerate(hmm_groups_used):
            for hmm_id in all_HMM_ids[hmm_group]:
                hmm_id_to_hmm_group_name[hmm_id] = hmm_group
                hmm_id_to_hmm_group_index[hmm_id] = hmm_group_i
        
        # Group loop
        hit_info_by_genome_feature_and_hmm = dict()  # used to build DomainAnnotation object
        total_hit_cnts = dict()
        accepted_hit_cnts = dict()
        hit_cnt_by_genome_and_model = dict()
        hit_genes_by_genome_and_model = dict()
        hit_accept_something = dict()
        output_hit_TAB_file_paths = dict()
        #output_hit_MSA_file_paths = dict()
        objects_created_refs_coalesce = dict()
        objects_created_refs_by_hmm_id = dict()

        for hmm_group_i, hmm_group in enumerate(all_HMM_groups_order):
            self.log(console, "PROCESSING HMM GROUP: " + hmm_group)  # DEBUG

            hit_accept_something[hmm_group] = False

            if hmm_group not in hmm_groups_used:
                for hmm_id in all_HMM_ids[hmm_group]:
                    total_hit_cnts[hmm_id] = 0
                    accepted_hit_cnts[hmm_id] = 0
                continue

            ## iterate through HMMs and scan input_many DBs
            #
            output_filtered_fasta_file_paths = []
            output_hits_flags = []
            coalesced_sequenceObjs = []
            coalesce_featureIds_element_ordering = []
            coalesce_featureIds_genome_ordering = []
            html_report_chunks = []

            # HMM loop
            for hmm_i, hmm_id in enumerate(input_HMM_ids[hmm_group]):

                self.log(console, "PROCESSING HMM: " + hmm_id)  # DEBUG

                # init hit counts
                total_hit_cnts[hmm_id] = 0
                accepted_hit_cnts[hmm_id] = 0
                html_report_chunks.append(None)

                # set paths
                #
                hmmer_dir = os.path.join(self.output_dir, hmm_id)  # this must match above
                if not os.path.exists(hmmer_dir):
                    os.makedirs(hmmer_dir)
                HMM_file_path = os.path.join(hmmer_dir, hmm_id + ".hmm")

                # create HMM file
                with open(HMM_file_path, 'w', 0) as hmm_handle:
                    hmm_handle.write("\n".join(HMM_bufs[hmm_id]) + "\n")

                if not os.path.isfile(HMM_file_path):
                    raise ValueError("HMMER_BUILD failed to create HMM file '" + HMM_file_path + "'")
                elif not os.path.getsize(HMM_file_path) > 0:
                    raise ValueError("HMMER_BUILD created empty HMM file '" + HMM_file_path + "'")

                # get model len (this can also be obtained from the -domtblout output)
                model_len[hmm_id] = 0
                for HMM_line in HMM_bufs[hmm_id]:
                    if HMM_line.startswith('LENG '):
                        model_len[hmm_id] = int(HMM_line.replace('LENG ','').strip())
                        break
                if model_len[hmm_id] == 0:
                    raise ValueError ("No length found in HMM file for model "+hmm_id)
                

                ### Construct the HMMER_SEARCH command for each target
                #
                # SYNTAX (from http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf)
                #
                # hmmsearch [--cut_tc] --domtblout <TAB_out> -A <MSA_out> --noali --notextw -E <e_value> -T <bit_score> <hmmfile> <seqdb>
                #
                for input_many_ref in input_many_refs:
                    hmmer_search_bin = self.HMMER_SEARCH
                    hmmer_search_cmd = [hmmer_search_bin]

                    # check for necessary files
                    if not os.path.isfile(hmmer_search_bin):
                        raise ValueError("no such file '" + hmmer_search_bin + "'")
                    if not os.path.isfile(HMM_file_path):
                        raise ValueError("no such file '" + HMM_file_path + "'")
                    elif not os.path.getsize(HMM_file_path):
                        raise ValueError("empty file '" + HMM_file_path + "'")
                    if not os.path.isfile(many_forward_reads_file_paths[input_many_ref]):
                        raise ValueError("no such file '" + many_forward_reads_file_paths[input_many_ref] + "'")
                    elif not os.path.getsize(many_forward_reads_file_paths[input_many_ref]):
                        raise ValueError("empty file '" + many_forward_reads_file_paths[input_many_ref] + "'")

                    output_hit_TAB_file_path = os.path.join(hmmer_dir, hmm_id +'-'+ input_many_names[input_many_ref] + '.hitout.txt')
                    #output_hit_MSA_file_path = os.path.join(hmmer_dir, hmm_id + '.msaout.txt')
                    #output_filtered_fasta_file_path = os.path.join(hmmer_dir, hmm_id + '.output_filtered.fasta')
                    if hmm_id not in output_hit_TAB_file_paths:
                        output_hit_TAB_file_paths[hmm_id] = dict()
                    output_hit_TAB_file_paths[hmm_id][input_many_ref] = output_hit_TAB_file_path
                    #output_hit_MSA_file_paths[hmm_id] = output_hit_MSA_file_path
                    #output_filtered_fasta_file_paths.append(output_filtered_fasta_file_path)

                    # this is command for basic search mode
                    hmmer_search_cmd.append('--domtblout')
                    hmmer_search_cmd.append(output_hit_TAB_file_path)
                    #hmmer_search_cmd.append('-A')
                    #hmmer_search_cmd.append(output_hit_MSA_file_path)
                    hmmer_search_cmd.append('--noali')
                    hmmer_search_cmd.append('--notextw')
                    if int(params.get('use_model_specific_thresholds',0)) == 1:
                        hmmer_search_cmd.append('--cut_tc')
                    else:
                        hmmer_search_cmd.append('--domE')  # can't use -T with -E, so we'll use -E
                        hmmer_search_cmd.append(str(params['e_value']))
                    hmmer_search_cmd.append(HMM_file_path)
                    hmmer_search_cmd.append(many_forward_reads_file_paths[input_many_ref])

                    # options
                    #if 'maxaccepts' in params:
                    #    if params['maxaccepts']:
                    #        hmmer_search_cmd.append('-max_target_seqs')
                    #        hmmer_search_cmd.append(str(params['maxaccepts']))
                    
                    # Run HMMER, capture output as it happens
                    #
                    #self.log(console, 'RUNNING HMMER_SEARCH:')
                    #self.log(console, '    '+' '.join(hmmer_search_cmd))
                    #report += "\n"+'running HMMER_SEARCH:'+"\n"
                    #report += '    '+' '.join(hmmer_search_cmd)+"\n"
                    
                    self.log (console, "RUNNING HMMER against target "+input_many_names[input_many_ref]) # DEBUG

                    # DEBUG
                    #self.log(console, "\n".join(hmmer_search_cmd))
                    
                    p = subprocess.Popen(hmmer_search_cmd,
                                         cwd=self.output_dir,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT,
                                         shell=False)
                    
                    while True:
                        line = p.stdout.readline()
                        if not line:
                            break
                        #self.log(console, line.replace('\n', ''))  # DEBUG

                    p.stdout.close()
                    p.wait()
                    #self.log(console, 'return code: ' + str(p.returncode))
                    if p.returncode != 0:
                        raise ValueError('Error running HMMER_SEARCH, return code: ' + str(p.returncode) +
                                         '\n\n' + '\n'.join(console))

                    #self.log (console, "RAN HMMER") # DEBUG

                    # Check for output
                    if not os.path.isfile(output_hit_TAB_file_path):
                        raise ValueError("HMMER_SEARCH failed to create TAB file '" + output_hit_TAB_file_path + "'")
                    elif not os.path.getsize(output_hit_TAB_file_path) > 0:
                        raise ValueError("HMMER_SEARCH created empty TAB file '" + output_hit_TAB_file_path + "'")
                    """
                    if not os.path.isfile(output_hit_MSA_file_path):
                    raise ValueError("HMMER_SEARCH failed to create MSA file '" + output_hit_MSA_file_path + "'")
                    elif not os.path.getsize(output_hit_MSA_file_path) > 0:
                    #raise ValueError("HMMER_SEARCH created empty MSA file '"+output_hit_MSA_file_path+"'")
                    #self.log(console,"HMMER_SEARCH created empty MSA file '"+output_hit_MSA_file_path+"'")
                    self.log(console, "\tHMMER_SEARCH: No hits")
                    continue
                    """
                
                    # DEBUG
                    #self.log(console, "DEBUG: output_hit_TAB_file_path: '"+str(output_hit_TAB_file_path))
                    #self.log(console, "DEBUG: output_hit_MSA_file_path: '"+str(output_hit_MSA_file_path))
                    # DEBUG
                    """
                    report = "MODEL ID:"+hmm_id+"\n"
                    report = "TAB:\n\n"
                    with open (output_hit_TAB_file_path, 'r') as output_handle:
                    for line in output_handle:
                    report += line+"\n"
                    report += "\n\nMSA:\n\n"
                    with open (output_hit_MSA_file_path, 'r') as output_handle:
                    for line in output_handle:
                    report += line+"\n"
                    self.log(console, report)
                    """

                ### Parse the HMMER tabular output and store ids to filter many set to make filtered object to save back to KBase
                #
                #self.log(console, 'PARSING HMMER SEARCH TAB OUTPUT')
                hits_found = False
                for input_many_ref in input_many_refs:
                    with open(output_hit_TAB_file_paths[hmm_id][input_many_ref], "r") as output_hit_TAB_file_handle:
                        output_aln_buf = output_hit_TAB_file_handle.readlines()

                    # check for hits
                    for line in output_aln_buf:
                        if line.startswith('#'):
                            continue
                        hits_found = True
                        break

                if not hits_found:
                    self.log(console, "\tHMMER_SEARCH: No hits")
                    continue
                
                self.log (console, "PARSING HMMER OUTPUT") # DEBUG
                    
                hit_seq_ids = dict()
                accepted_hit_cnt = 0  # may be more than one hit per gene
                accepted_hit_seq_ids = dict()
                accept_fids = dict()
                filtering_fields = dict()
                hit_order = []
                hit_buf = dict()
                hit_source_ref  = dict()
                hit_source_type = dict()
                
                for input_many_ref in input_many_refs:
                    many_type_name = many_type_names[input_many_ref]
                    with open(output_hit_TAB_file_paths[hmm_id][input_many_ref], "r") as output_hit_TAB_file_handle:
                        output_aln_buf = output_hit_TAB_file_handle.readlines()

                    for line in output_aln_buf:
                        if line.startswith('#'):
                            #if not header_done:
                            #    hit_buf.append(line)
                            continue
                        #header_done = True
                        #self.log(console,'HIT LINE: '+line)  # DEBUG
                        """ format for -tblout.  we're now using -domtblout" 
                        hit_info = re.split('\s+', line)
                        hit_seq_id = hit_info[0]
                        hit_accession = hit_info[1]
                        query_name = hit_info[2]
                        query_accession = hit_info[3]
                        hit_e_value = float(hit_info[4])
                        hit_bitscore = float(hit_info[5])
                        hit_bias = float(hit_info[6])
                        hit_e_value_best_dom = float(hit_info[7])
                        hit_bitscore_best_dom = float(hit_info[8])
                        hit_bias_best_dom = float(hit_info[9])
                        hit_expected_dom_n = float(hit_info[10])
                        hit_regions = float(hit_info[11])
                        hit_regions_multidom = float(hit_info[12])
                        hit_overlaps = float(hit_info[13])
                        hit_envelopes = float(hit_info[14])
                        hit_dom_n = float(hit_info[15])
                        hit_doms_within_rep_thresh = float(hit_info[16])
                        hit_doms_within_inc_thresh = float(hit_info[17])
                        hit_desc = hit_info[18]
                        """
                        # format for -domtblout http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf
                        hit_info = re.split('\s+', line)
                        hit_seq_id = hit_info[0]
                        hit_accession = hit_info[1]
                        hit_seq_len = int(hit_info[2])  # diff from -tblout format
                        query_name = hit_info[3]
                        query_accession = hit_info[4]
                        query_model_len = int(hit_info[5])  # diff from -tblout format
                        hit_e_value_alldoms = float(hit_info[6])
                        hit_bitscore_alldoms = float(hit_info[7])
                        hit_bias_alldoms = float(hit_info[8])
                        hit_dom_n = int(hit_info[9])  # diff from -tblout format
                        hit_dom_total = int(hit_info[10])  # diff from -tblout format
                        hit_c_e_value_this_dom = float(hit_info[11])  # diff from -tblout format
                        hit_i_e_value_this_dom = float(hit_info[12])  # diff from -tblout format
                        hit_bitscore_this_dom = float(hit_info[13])  # diff from -tblout format
                        hit_bias_this_dom = float(hit_info[14])  # diff from -tblout format
                        query_beg_pos = int(hit_info[15])  # diff from -tblout format
                        query_end_pos = int(hit_info[16])  # diff from -tblout format
                        hit_beg_pos = int(hit_info[17])  # diff from -tblout format
                        hit_end_pos = int(hit_info[18])  # diff from -tblout format
                        envelope_beg_pos = int(hit_info[19])  # diff from -tblout format
                        envelope_end_pos = int(hit_info[20])  # diff from -tblout format
                        hit_acc_posterior_prob = float(hit_info[21])  # diff from -tblout format
                        hit_desc = hit_info[22]

                        if hit_seq_id not in hit_seq_ids:
                            hit_seq_ids[hit_seq_id] = True
                            hit_order.append(hit_seq_id)
                            hit_buf[hit_seq_id] = []
                            hit_source_ref[hit_seq_id] = input_many_ref
                            hit_source_type[hit_seq_id] = many_type_name
                            filtering_fields[hit_seq_id] = []

                        # store all non-filtered domain hits for domainAnnotation
                        #hit_info_by_genome_feature_and_hmm = dict()
                        [hit_genome_ref, hit_feature_id] = self._parse_genome_and_feature_id_from_hit_id (hit_seq_id, many_type_names[input_many_ref], input_many_ref, genome_id_feature_id_delim)
                        hit_alnlen = query_end_pos - query_beg_pos + 1

                        if hit_genome_ref not in hit_info_by_genome_feature_and_hmm:
                            hit_info_by_genome_feature_and_hmm[hit_genome_ref] = dict()
                        if hit_feature_id not in hit_info_by_genome_feature_and_hmm[hit_genome_ref]:
                            hit_info_by_genome_feature_and_hmm[hit_genome_ref][hit_feature_id] = dict()
                        if hmm_id not in hit_info_by_genome_feature_and_hmm[hit_genome_ref][hit_feature_id]:
                            hit_info_by_genome_feature_and_hmm[hit_genome_ref][hit_feature_id][hmm_id] = list()

                        filter = False
                        this_filter_fields = dict()
                        #self.log(console,"HIT_SEQ_ID: '"+hit_seq_id+"'")
                        #if 'ident_thresh' in params and float(params['ident_thresh']) > float(high_bitscore_ident[hit_seq_id]):
                        #    continue
                        if 'bitscore' in params and float(params['bitscore']) > float(hit_bitscore_this_dom):
                            self.log(console, "model "+hmm_id+" filtering "+hit_seq_id+" with bitscore "+str(hit_bitscore_this_dom))  # DEBUG
                            filter = True
                            this_filter_fields['bitscore'] = True
                        if 'model_cov_perc' in params and float(params['model_cov_perc']) > 100.0 * float(hit_alnlen) / float(model_len[hmm_id]):
                            self.log(console, "model "+hmm_id+" filtering "+hit_seq_id+" with model_cov "+str(hit_alnlen))  # DEBUG
                            filter = True
                            this_filter_fields['model_cov_perc'] = True
                        if 'maxaccepts' in params and params['maxaccepts'] != None and accepted_hit_cnt >= int(params['maxaccepts']):
                            self.log(console, "model "+hmm_id+" filtering "+hit_seq_id+" with maxaccepts "+str(accepted_hit_cnt))  # DEBUG
                            filter = True
                            this_filter_fields['maxaccepts'] = True

                        filtering_fields[hit_seq_id].append(this_filter_fields)
                        if filter:
                            continue

                        # store hit
                        accepted_hit_cnt += 1
                        accepted_hit_seq_ids[hit_seq_id] = True
                        hit_buf[hit_seq_id].append(line)
                        hit_accept_something[hmm_group] = True
                        #self.log(console, "HIT: '"+hit_seq_id+"'")  # DEBUG

                        # capture accepted hit count by genome_ref and model
                        if many_type_names[input_many_ref] == 'Genome' \
                           or many_type_names[input_many_ref] == 'AnnotatedMetagenomeAssembly':
                            genome_ref = input_many_ref
                            gene_id = hit_seq_id
                        else:
                            [genome_ref, gene_id] = hit_seq_id.split(genome_id_feature_id_delim)
                        #self.log(console, "DEBUG: genome_ref: '"+str(genome_ref)+"'")
                        #self.log(console, "DEBUG: input_hmm_name: '"+str(hmm_id)+"'")
                        if genome_ref not in hit_cnt_by_genome_and_model:
                            hit_cnt_by_genome_and_model[genome_ref] = dict()
                            hit_genes_by_genome_and_model[genome_ref] = dict()
                        if hmm_id not in hit_cnt_by_genome_and_model[genome_ref]:
                            hit_cnt_by_genome_and_model[genome_ref][hmm_id] = 0
                            hit_genes_by_genome_and_model[genome_ref][hmm_id] = []
                        hit_cnt_by_genome_and_model[genome_ref][hmm_id] += 1
                        hit_genes_by_genome_and_model[genome_ref][hmm_id].append(gene_id)

                    
                        # prep for storing DomainAnnotation
                        # domain_place: tuple<int start_in_feature,int stop_in_feature,float evalue, float bitscore,float domain_coverage>).
                        this_hit = {'start_in_feature': hit_beg_pos,
                                    'stop_in_feature': hit_end_pos,
                                    'evalue': hit_i_e_value_this_dom,
                                    'bitscore': hit_bitscore_this_dom,
                                    'domain_coverage': round (hit_alnlen / model_len[hmm_id],2)
                                
                        }
                        hit_info_by_genome_feature_and_hmm[hit_genome_ref][hit_feature_id][hmm_id].append(this_hit)
                        
                    
                ### Create output objects
                #
                total_hit_cnts[hmm_id] = len(hit_order)  # this is just gene total
                accepted_hit_cnts[hmm_id] = accepted_hit_cnt  # this may incude multiple hits per gene

                if accepted_hit_cnt == 0:
                    self.log(console, "\tNO ACCEPTED HITS ABOVE FILTERS")
                elif (not params.get('save_ALL_featureSets') or int(params['save_ALL_featureSets']) != 1) and \
                     not explicitly_requested_models.get(hmm_id):
                    self.log(console, "\tMODEL "+hmm_id+" NOT EXPLICITLY REQUESTED, SO NOT SAVING FEATURESET.  MUST EXPLICITLY REQUEST MODEL OR CHANGE 'save_ALL_featureSets' to TRUE.")
                elif int(params.get('save_ANY_featureSets',1)) != 1:
                    self.log(console, "\tsave_ANY_featureSets set to FALSE. SO NOT SAVING FEATURESET.  MUST CHANGE 'save_ANY_featureSets' to TRUE.")
                else:
                    self.log(console, "\tEXTRACTING ACCEPTED HITS FROM INPUT for model "+hmm_id)
                    ##self.log(console, 'MANY_TYPE_NAME: '+many_type_name)  # DEBUG
                    output_featureSet = dict()
                    output_featureSet['description'] = search_tool_name + "_Search filtered"
                    output_featureSet['element_ordering'] = []
                    output_featureSet['elements'] = dict()

                    # go through input targets
                    for input_many_ref in input_many_refs:

                        many_type_name = many_type_names[input_many_ref]
                        
                        # SequenceSet input -> SequenceSet output
                        #
                        if many_type_name == 'SequenceSet':
                            output_sequenceSet = dict()

                            if 'sequence_set_id' in input_many_sequenceSet and input_many_sequenceSet['sequence_set_id'] != None:
                                output_sequenceSet['sequence_set_id'] = input_many_sequenceSet['sequence_set_id'] + \
                                                                        "." + search_tool_name + "_Search_filtered"
                            else:
                                output_sequenceSet['sequence_set_id'] = search_tool_name + "_Search_filtered"
                            if 'description' in input_many_sequenceSet and input_many_sequenceSet['description'] != None:
                                output_sequenceSet['description'] = input_many_sequenceSet['description'] + \
                                                                    " - " + search_tool_name + "_Search filtered"
                            else:
                                output_sequenceSet['description'] = search_tool_anme + "_Search filtered"

                            #self.log(console,"ADDING SEQUENCES TO SEQUENCESET")
                            output_sequenceSet['sequences'] = []

                            for seq_obj in input_many_sequenceSet['sequences']:
                                header_id = seq_obj['sequence_id']
                                #header_desc = seq_obj['description']
                                #sequence_str = seq_obj['sequence']

                                id_untrans = header_id
                                id_trans = re.sub('\|', ':', id_untrans)
                                if id_trans in accepted_hit_seq_ids or id_untrans in accepted_hit_seq_ids:
                                    #self.log(console, 'FOUND HIT '+header_id)  # DEBUG
                                    accept_fids[id_untrans] = True
                                    output_sequenceSet['sequences'].append(seq_obj)

                        # FeatureSet input -> FeatureSet output
                        #
                        elif many_type_name == 'FeatureSet':
                            input_many_featureSet = input_many_objs[input_many_ref]
                            fId_list = input_many_featureSet['elements'].keys()
                            #self.log(console,"ADDING FEATURES TO FEATURESET")
                            for fId in sorted(fId_list):
                                for genome_ref in input_many_featureSet['elements'][fId]:
                                    id_untrans = genome_ref + genome_id_feature_id_delim + fId
                                    id_trans = re.sub('\|', ':', id_untrans)
                                    if id_trans in accepted_hit_seq_ids or id_untrans in accepted_hit_seq_ids:
                                        #self.log(console, 'FOUND HIT '+fId)  # DEBUG
                                        accept_fids[id_untrans] = True
                                        #fId = id_untrans  # don't change fId for output FeatureSet
                                        if fid not in output_featureSet['elements']:
                                            output_featureSet['elements'][fId] = []
                                            output_featureSet['element_ordering'].append(fId)
                                        if genome_ref not in output_featureSet['elements'][fId]:
                                            output_featureSet['elements'][fId].append(genome_ref)

                        # Parse Genome hits into FeatureSet
                        #
                        elif many_type_name == 'Genome':
                            genome_ref = input_many_ref
                            for fid in feature_ids[input_many_ref]['basic']:
                                id_untrans = fid
                                id_trans = re.sub('\|', ':', id_untrans)
                                if id_trans in accepted_hit_seq_ids or id_untrans in accepted_hit_seq_ids:
                                    #self.log(console, 'FOUND HIT '+fid)  # DEBUG
                                    #output_featureSet['element_ordering'].append(fid)
                                    accept_fids[id_untrans] = True
                                    #fid = input_many_ref+genome_id_feature_id_delim+id_untrans  # don't change fId for output FeatureSet
                                    if fid not in output_featureSet['elements']:
                                        output_featureSet['element_ordering'].append(fid)
                                        output_featureSet['elements'][fid] = []
                                    if genome_ref not in output_featureSet['elements'][fid]:
                                        output_featureSet['elements'][fid].append(genome_ref)

                        # Parse GenomeSet hits into FeatureSet
                        #
                        elif many_type_name == 'GenomeSet':
                            #self.log(console,"READING HITS FOR GENOMES")  # DEBUG
                            input_many_genomeSet = input_many_objs[input_many_ref]
                            for genome_id in feature_ids[input_many_ref]['by_genome_id'].keys():
                                #self.log(console,"READING HITS FOR GENOME "+genome_id)  # DEBUG
                                genome_ref = input_many_genomeSet['elements'][genome_id]['ref']
                                for feature_id in feature_ids[input_many_ref]['by_genome_id'][genome_id]:
                                    id_untrans = genome_ref + genome_id_feature_id_delim + feature_id
                                    id_trans = re.sub('\|', ':', id_untrans)
                                    if id_trans in accepted_hit_seq_ids or id_untrans in accepted_hit_seq_ids:
                                        #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                                        #output_featureSet['element_ordering'].append(feature['id'])
                                        accept_fids[id_untrans] = True
                                        #feature_id = id_untrans  # don't change fId for output FeatureSet
                                        if feature_id not in output_featureSet['elements']:
                                            output_featureSet['elements'][feature_id] = []
                                            output_featureSet['element_ordering'].append(feature_id)
                                        if genome_ref not in output_featureSet['elements'][feature_id]:
                                            output_featureSet['elements'][feature_id].append(genome_ref)

                        # Parse SpeciesTree hits into FeatureSet
                        #
                        elif many_type_name == 'Tree':
                            #self.log(console,"READING HITS FOR GENOMES")  # DEBUG
                            input_many_tree = input_many_objs[input_many_ref]
                            for genome_id in feature_ids[input_many_ref]['by_genome_id'].keys():
                                #self.log(console,"READING HITS FOR GENOME "+genome_id)  # DEBUG
                                genome_ref = input_many_tree['ws_refs'][genome_id]['g'][0]
                                for feature_id in feature_ids[input_many_ref]['by_genome_id'][genome_id]:
                                    id_untrans = genome_ref + genome_id_feature_id_delim + feature_id
                                    id_trans = re.sub('\|', ':', id_untrans)
                                    if id_trans in accepted_hit_seq_ids or id_untrans in accepted_hit_seq_ids:
                                        #self.log(console, 'FOUND HIT: '+feature['id'])  # DEBUG
                                        #output_featureSet['element_ordering'].append(feature['id'])
                                        accept_fids[id_untrans] = True
                                        #feature_id = id_untrans  # don't change fId for output FeatureSet
                                        if feature_id not in output_featureSet['elements']:
                                            output_featureSet['elements'][feature_id] = []
                                            output_featureSet['element_ordering'].append(feature_id)
                                        if genome_ref not in output_featureSet['elements'][feature_id]:
                                            output_featureSet['elements'][feature_id].append(genome_ref)

                        # Parse AnnotatedMetagenomeAssembly hits into FeatureSet
                        #
                        elif many_type_name == 'AnnotatedMetagenomeAssembly':
                            #seq_total = 0
                            for fid in feature_ids[input_many_ref]['basic']:
                                #if fid == 'AWN69_RS07145' or fid == 'AWN69_RS13375':
                                #    self.log(console, 'CHECKING FID '+fid)  # DEBUG
                                #seq_total += 1
                                id_untrans = fid
                                id_trans = re.sub ('\|',':',id_untrans)
                                #print ("TESTING FEATURES: ID_UNTRANS: '"+id_untrans+"'")  # DEBUG
                                #print ("TESTING FEATURES: ID_TRANS: '"+id_trans+"'")  # DEBUG
                                if id_trans in accepted_hit_seq_ids or id_untrans in accepted_hit_seq_ids:
                                    self.log(console, 'FOUND HIT '+fid)  # DEBUG
                                    #output_featureSet['element_ordering'].append(fid)
                                    accept_fids[id_untrans] = True
                                    #fid = input_many_ref+self.genome_id_feature_id_delim+id_untrans  # don't change fId for output FeatureSet
                                    ama_ref = input_many_ref
                                    if fid not in output_featureSet['elements']:
                                        output_featureSet['element_ordering'].append(fid)
                                        output_featureSet['elements'][fid] = []
                                    if ama_ref not in output_featureSet['elements'][fid]:
                                        output_featureSet['elements'][fid] = [ama_ref]

                                    
                    ### Create output object
                    #
                    if 'coalesce_output' in params and int(params['coalesce_output']) == 1:
                        if len(invalid_msgs) == 0:
                            if len(accepted_hit_seq_ids.keys()) == 0:   # Note, this is after filtering, so there may be more unfiltered hits
                                #self.log(console,"No Object to Upload for HMM "+hmm_id)  # DEBUG
                                continue

                            # accumulate hits into coalesce object
                            #
                            #if many_type_name == 'SequenceSet':  # input many SequenceSet -> save SequenceSet
                            #    for seq_obj in output_sequenceSet['sequences']:
                            #        coalesced_sequenceObjs.append(seq_obj)
                            #
                            #else:  # input FeatureSet, Genome, and GenomeSet -> upload FeatureSet output
                            for fId in output_featureSet['element_ordering']:
                                coalesce_featureIds_element_ordering.append(fId)
                                #coalesce_featureIds_genome_ordering.append(output_featureSet['elements'][fId][0])
                                for this_genome_ref in output_featureSet['elements'][fId]:
                                    coalesce_featureIds_genome_ordering.append(this_genome_ref)

                    else:  # keep output separate  Upload results if coalesce_output is 0
                        hmm_obj_id = hmm_id
                        if '+' in hmm_id:
                            hmm_obj_id = hmm_obj_id.replace('+','-')
                        output_name = hmm_obj_id + '-' + params['output_filtered_name']

                        if len(invalid_msgs) == 0:
                            if len(accepted_hit_seq_ids.keys()) == 0:   # Note, this is after filtering, so there may be more unfiltered hits
                                #self.log(console,"No Object to Upload for HMM "+hmm_id)  # DEBUG
                                continue

                            #self.log(console,"Uploading results Object HMM "+hmm_id)  # DEBUG

                            # input many SequenceSet -> save SequenceSet
                            #
                            '''
                            if many_type_name == 'SequenceSet':
                                new_obj_info = self.wsClient.save_objects({
                                    'workspace': params['workspace_name'],
                                    'objects': [{
                                        'type': 'KBaseSequences.SequenceSet',
                                        'data': output_sequenceSet,
                                        'name': output_name,
                                        'meta': {},
                                        'provenance': provenance
                                    }]
                                })[0]
                            '''
                            #else:  # input FeatureSet, Genome, GenomeSet, Tree, and AMA -> upload FeatureSet output
                            new_obj_info = self.wsClient.save_objects({
                                'workspace': params['workspace_name'],
                                'objects': [{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_featureSet,
                                    'name': output_name,
                                    'meta': {},
                                    'provenance': provenance
                                }]
                            })[0]

                            [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I,
                                WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
                            objects_created_refs_by_hmm_id[hmm_id] = str(
                                new_obj_info[WSID_I]) + '/' + str(new_obj_info[OBJID_I])

                #### Build output report chunks
                ##
                #self.log(console,"BUILDING REPORT CHUNK for HMM["+str(hmm_i)+"] "+hmm_id)  # DEBUG
                if len(invalid_msgs) == 0:

                    # text report
                    #
                    report += 'HMM[' + str(hmm_i) + ']: ' + hmm_id + "\n"
                    report += 'sequences in search db: ' + str(seq_total) + "\n"
                    report += 'sequences in hit set: ' + str(total_hit_cnts[hmm_id]) + "\n"
                    report += 'domain hits in accepted hit set: ' + str(accepted_hit_cnts[hmm_id]) + "\n"
                    report += "\n"
                    #for line in hit_buf:
                    #    report += line
                    #self.log (console, report)

                    # build html report chunk
                    feature_id_to_function = dict()
                    genome_ref_to_obj_name = dict()
                    ama_ref_to_obj_name    = dict()
                    genome_ref_to_sci_name = dict()
                    for input_many_ref in input_many_refs:
                        many_type_name = many_type_names[input_many_ref]

                        print ("SETTING FEATURE_ID_TO_FUNCTION")  # DEBUG
                        for genome_ref in obj2file_retVal[input_many_ref]['feature_id_to_function'].keys():
                            if genome_ref not in feature_id_to_function:
                                feature_id_to_function[genome_ref] = dict()
                            for fid in  obj2file_retVal[input_many_ref]['feature_id_to_function'][genome_ref].keys():
                                feature_id_to_function[genome_ref][fid] = obj2file_retVal[input_many_ref]['feature_id_to_function'][genome_ref][fid]
                        if many_type_name == 'Genome' or many_type_name == 'GenomeSet' or many_type_name == 'FeatureSet' or many_type_name == 'Tree':
                            for genome_ref in obj2file_retVal[input_many_ref]['genome_ref_to_obj_name'].keys():
                                genome_ref_to_obj_name[genome_ref] = obj2file_retVal[input_many_ref]['genome_ref_to_obj_name'][genome_ref]
                                genome_ref_to_sci_name[genome_ref] = obj2file_retVal[input_many_ref]['genome_ref_to_sci_name'][genome_ref]
                        elif many_type_name == 'AnnotatedMetagenomeAssembly':
                            for ama_ref in obj2file_retVal[input_many_ref]['ama_ref_to_obj_name'].keys():
                                #self.log(console, "AMA_REF: "+ama_ref)
                                #ama_ref_to_obj_name[ama_ref] = AnnotatedMetagenomeAssemblyToFASTA_retVal['ama_ref_to_obj_name'][ama_ref]
                                ama_ref_to_obj_name[ama_ref] = obj2file_retVal[input_many_ref]['ama_ref_to_obj_name'][ama_ref]

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

                    for hit_seq_id in hit_order:

                        input_many_ref = hit_source_ref[hit_seq_id]
                        many_type_name = hit_source_type[hit_seq_id]

                        for hit_i,line in enumerate(hit_buf[hit_seq_id]):
                            line = line.strip()
                            if line == '' or line.startswith('#'):
                                continue

                            """ This is -tbleout format
                            [hit_id, hit_accession, query_name, query_accession, e_value, bit_score, bias, e_value_best_dom, bit_score_best_dom, bias_best_dom, expected_dom_n,
                            regions, regions_multidom, overlaps, envelopes, dom_n, doms_within_rep_thresh, doms_within_inc_thresh, hit_desc] = re.split('\s+', line)[0:19]
                            """

                            # use -domtbleout format
                            #                [query_id, hit_id, identity, aln_len, mismatches, gap_openings, q_beg, q_end, h_beg, h_end, e_value, bit_score] = line.split("\t")[0:12]
                            #                identity = str(round(float(identity), 1))
                            #                if identity == '100.0':  identity = '100'
                            [hit_seq_id,
                             hit_accession,
                             hit_seq_len,
                             query_name,
                             query_accession,
                             query_model_len,
                             hit_e_value_alldoms,
                             hit_bitscore_alldoms,
                             hit_bias_alldoms,
                             hit_dom_n,
                             hit_dom_total,
                             hit_c_e_value_this_dom,
                             hit_i_e_value_this_dom,
                             hit_bitscore_this_dom,
                             hit_bias_this_dom,
                             query_beg_pos,
                             query_end_pos,
                             hit_beg_pos,
                             hit_end_pos,
                             envelope_beg_pos,
                             envelope_end_pos,
                             hit_acc_posterior_prob,
                             hit_desc
                            ] = re.split('\s+', line)[0:23]
                        
                            # get coords with respect to hit sequence and model
                            h_len = int(hit_seq_len)
                            h_beg = int(hit_beg_pos)
                            h_end = int(hit_end_pos)
                            e_beg = int(envelope_beg_pos)
                            e_end = int(envelope_end_pos)
                            q_len = int(query_model_len)
                            q_beg = int(query_beg_pos)
                            q_end = int(query_end_pos)
                            aln_len = abs(q_end - q_beg) + 1
                            aln_len_perc = round(100.0 * float(aln_len) / float(q_len), 1)


                            #if many_type_name == 'SingleEndLibrary':
                            #    pass
                            #elif many_type_name == 'SequenceSet':
                            if many_type_name == 'SequenceSet':
                                pass
                            elif many_type_name == 'Genome' or \
                                 many_type_name == 'AnnotatedMetagenomeAssembly' or \
                                 many_type_name == 'GenomeSet' or \
                                 many_type_name == 'Tree' or \
                                 many_type_name == 'FeatureSet':
                                
                                if 'Set' in many_type_name or many_type_name == 'Tree':
                                    [genome_ref, hit_fid] = hit_seq_id.split(genome_id_feature_id_delim)
                                else:
                                    genome_ref = input_many_ref
                                    hit_fid = hit_seq_id

                                # can't just use hit_fid because may have pipes translated and can't translate back
                                fid_lookup = None
                                for fid in feature_id_to_function[genome_ref].keys():
                                    id_untrans = fid
                                    id_trans = re.sub('\|', ':', id_untrans)

                                    #self.log (console, "SCANNING FIDS.  HIT_FID: '"+str(hit_fid)+"' FID: '"+str(fid)+"' TRANS: '"+str(id_trans)+"'")  # DEBUG

                                    if id_untrans == hit_fid or id_trans == hit_fid:
                                        #self.log (console, "GOT ONE!")  # DEBUG
                                        if many_type_name == 'Genome' or many_type_name == 'AnnotatedMetagenomeAssembly':
                                            accept_id = fid
                                        elif many_type_name == 'GenomeSet' or many_type_name == 'FeatureSet' or many_type_name == 'Tree':
                                            accept_id = genome_ref + genome_id_feature_id_delim + fid
                                        if accept_id in accept_fids:
                                            row_color = accept_row_color
                                        else:
                                            row_color = reject_row_color
                                        fid_lookup = fid
                                        break
                                    #self.log (console, "HIT_FID: '"+str(hit_fid)+"' FID_LOOKUP: '"+str(fid_lookup)+"'")  # DEBUG
                                if fid_lookup == None:
                                    raise ValueError("unable to find fid for hit_fid: '" + str(hit_fid))
                                elif fid_lookup not in feature_id_to_function[genome_ref]:
                                    raise ValueError("unable to find function for fid: '" + str(fid_lookup))
                                fid_disp = re.sub(r"^.*\.([^\.]+)\.([^\.]+)$", r"\1.\2", fid_lookup)

                                func_disp = feature_id_to_function[genome_ref][fid_lookup]

                                # set genome_disp_name
                                genome_disp_name = self._get_genome_disp_name (params,
                                                                               genome_ref,
                                                                               many_type_name,
                                                                               ama_ref_to_obj_name,
                                                                               genome_ref_to_obj_name,
                                                                               genome_ref_to_sci_name)

                                # build html report table line
                                html_report_chunk += ['<tr bgcolor="' + row_color + '">']
                                #html_report_chunk += ['<tr bgcolor="'+'white'+'">']  # DEBUG
                                # add overlap bar

                                # coverage graphic (with respect to hit seq)
                                html_report_chunk += ['<td valign=middle align=center style="border-right:solid 1px ' +
                                                      border_body_color + '; border-bottom:solid 1px ' + border_body_color + '">']
                                html_report_chunk += ['<table style="height:' + str(bar_height) + 'px; width:' + str(
                                    bar_width) + 'px" border=0 cellpadding=0 cellspacing=0>']
                                full_len_pos = bar_width
                                aln_beg_pos = int(float(bar_width) * float(int(h_beg) - 1) / float(int(h_len) - 1))
                                aln_end_pos = int(float(bar_width) * float(int(h_end) - 1) / float(int(h_len) - 1))
                                cell_pix_height = str(int(round(float(bar_height) / 3.0, 0)))

                                cell_color = ['', '', '']
                                cell_width = []
                                cell_width.append(aln_beg_pos)
                                cell_width.append(aln_end_pos - aln_beg_pos)
                                cell_width.append(bar_width - aln_end_pos)

                                for row_i in range(3):
                                    html_report_chunk += ['<tr style="height:' + cell_pix_height + 'px">']
                                    unalign_color = row_color
                                    if row_i == 1:
                                        unalign_color = bar_line_color
                                    cell_color[0] = unalign_color
                                    cell_color[1] = bar_color
                                    cell_color[2] = unalign_color

                                    for col_i in range(3):
                                        cell_pix_width = str(cell_width[col_i])
                                        cell_pix_color = cell_color[col_i]
                                        html_report_chunk += ['<td style="height:' + cell_pix_height +
                                                              'px; width:' + cell_pix_width + 'px" bgcolor="' + cell_pix_color + '"></td>']
                                    html_report_chunk += ['</tr>']
                                html_report_chunk += ['</table>']
                                html_report_chunk += ['</td>']

                                # add other cells
                                # fid
                                html_report_chunk += ['<td style="border-right:solid 1px ' + border_body_color + '; border-bottom:solid 1px ' +
                                                      border_body_color + '"><font color="' + text_color + '" size=' + text_fontsize + '>' + str(fid_disp) + '</font></td>']
                                #                    html_report_chunk += ['<td style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(hit_accession)+'</font></td>']
                                # func
                                html_report_chunk += ['<td style="border-right:solid 1px ' + border_body_color + '; border-bottom:solid 1px ' +
                                                      border_body_color + '"><font color="' + text_color + '" size=' + text_fontsize + '>' + func_disp + '</font></td>']
                                # genome name
                                html_report_chunk += ['<td style="border-right:solid 1px ' + border_body_color + '; border-bottom:solid 1px ' +
                                                      border_body_color + '"><font color="' + text_color + '" size=' + text_fontsize + '>' + genome_disp_name + '</font></td>']
                                # ident
                                #                    if 'ident_thresh' in filtering_fields[hit_id]:
                                #                       this_cell_color = reject_cell_color
                                #                   else:
                                #                       this_cell_color = row_color
                                #                   html_report_chunk += ['<td align=center bgcolor="'+this_cell_color+'" style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+str(identity)+'%</font></td>']

                                # aln len
                                if 'model_cov_perc' in filtering_fields[hit_seq_id]:
                                    this_cell_color = reject_cell_color
                                else:
                                    this_cell_color = row_color
                                html_report_chunk += ['<td align=center bgcolor="' + str(this_cell_color) + '" style="border-right:solid 1px ' + border_body_color + '; border-bottom:solid 1px ' +
                                                      border_body_color + '"><font color="' + text_color + '" size=' + text_fontsize + '>' + str(aln_len) + ' (' + str(aln_len_perc) + '%)</font></td>']

                                # evalue
                                html_report_chunk += ['<td align=center style="border-right:solid 1px ' + border_body_color + '; border-bottom:solid 1px ' +
                                                      border_body_color + '"><font color="' + text_color + '" size=' + text_fontsize + '><nobr>' + str(hit_i_e_value_this_dom) + '</nobr></font></td>']

                                # bit score
                                if 'bitscore' in filtering_fields[hit_seq_id]:
                                    this_cell_color = reject_cell_color
                                else:
                                    this_cell_color = row_color
                                html_report_chunk += ['<td align=center bgcolor="' + str(this_cell_color) + '" style="border-right:solid 1px ' + border_body_color + '; border-bottom:solid 1px ' +
                                                      border_body_color + '"><font color="' + text_color + '" size=' + text_fontsize + '><nobr>' + str(hit_bitscore_this_dom) + '</nobr></font></td>']
                                # bias
                                #                    html_report_chunk += ['<td align=center style="border-right:solid 1px '+border_body_color+'; border-bottom:solid 1px '+border_body_color+'"><font color="'+text_color+'" size='+text_fontsize+'><nobr>'+str(bias)+'</nobr><br><nobr>('+str(bias_best_dom)+')</nobr></font></td>']

                                # aln coords only for hit seq
                                html_report_chunk += ['<td align=center style="border-right:solid 1px ' + border_body_color + '; border-bottom:solid 1px ' + border_body_color +
                                                      '"><font color="' + text_color + '" size=' + text_fontsize + '><nobr>' + str(h_beg) + '-' + str(h_end) + '</nobr></font></td>']

                                # close chunk
                                html_report_chunk += ['</tr>']

                                
                    # attach chunk
                    if hmm_id not in total_hit_cnts:
                        self.log(console, "MODEL "+hmm_id+" NOT REQUESTED.  NOT ADDING TO HTML HIT REPORT.")
                        html_report_chunk_str = '<tr><td colspan=table_col_width><blockquote><i>Model '+hmm_id+' not requested</i></td></tr>'
                    elif total_hit_cnts[hmm_id] == 0:
                        #self.log(console, "NO HITS FOR HMM["+str(hmm_i)+"] "+hmm_id+".  NOT ADDING TO HTML HIT REPORT.")
                        html_report_chunk_str = '<tr><td colspan=table_col_width><blockquote><i>no hits found</i></td></tr>'
                    else:
                        html_report_chunk_str = "\n".join(html_report_chunk)
                    html_report_chunks[hmm_i] = html_report_chunk_str
                    #self.log(console, "HTML_REPORT_CHUNK: '"+str(html_report_chunk_str)+"'")  # DEBUG


            #### Create and Upload output object for hmm_group if coalesce_output is true
            ##
            if 'coalesce_output' in params and int(params['coalesce_output']) == 1:
                output_name = hmm_group + '-' + params['output_filtered_name']

                if len(invalid_msgs) == 0:
                    if not hit_accept_something[hmm_group]:
                        self.log(console, "No Coalesced Hits Object to Upload for all HMMs in Group " + hmm_group)  # DEBUG

                    elif int(params.get('save_ANY_featureSets',1)) != 1:
                        self.log(console, "save_ANY_featureSets not set to TRUE")  # DEBUG

                    else:
                        self.log(console, "Uploading Coalesced Hits Object for HMM Group " + hmm_group)  # DEBUG

                        if many_type_name == 'SequenceSet':  # input many SequenceSet -> save SequenceSet

                            output_sequenceSet['sequences'] = coalesced_sequenceObjs
                            new_obj_info = self.wsClient.save_objects({
                                'workspace': params['workspace_name'],
                                'objects': [{
                                    'type': 'KBaseSequences.SequenceSet',
                                    'data': output_sequenceSet,
                                    'name': output_name,
                                    'meta': {},
                                    'provenance': provenance
                                }]
                            })[0]

                        else:  # input FeatureSet, Genome, GenomeSet, Tree, and AMA -> upload FeatureSet output

                            output_featureSet['element_ordering'] = coalesce_featureIds_element_ordering
                            output_featureSet['elements'] = dict()
                            for f_i, fId in enumerate(output_featureSet['element_ordering']):
                                output_featureSet['elements'][fId] = []
                                output_featureSet['elements'][fId].append(coalesce_featureIds_genome_ordering[f_i])

                            new_obj_info = self.wsClient.save_objects({
                                'workspace': params['workspace_name'],
                                'objects': [{
                                    'type': 'KBaseCollections.FeatureSet',
                                    'data': output_featureSet,
                                    'name': output_name,
                                    'meta': {},
                                    'provenance': provenance
                                }]
                            })[0]

                        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I,
                            WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
                        objects_created_refs_coalesce[hmm_group] = str(
                            new_obj_info[WSID_I]) + '/' + str(new_obj_info[OBJID_I])

                        
            #### Set paths for output HTML
            ##
            html_output_dir = os.path.join(self.output_dir, 'html_output')
            if not os.path.exists(html_output_dir):
                os.makedirs(html_output_dir)
            html_search_file = search_tool_name + '_Search-' + str(hmm_group_i) + '-' + str(hmm_group) + '.html'
            html_search_path = os.path.join(html_output_dir, html_search_file)
            html_profile_file = search_tool_name + '_Profile.html'
            html_profile_path = os.path.join(html_output_dir, html_profile_file)

            #### Build Search output report (and assemble html chunks)
            ##
            self.log(console, "BUILDING SEARCH REPORT ")  # DEBUG
            if len(invalid_msgs) == 0:

                # build html report
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
                html_report_lines += ['<title>KBase '+model_group_config['group_name'] + ' Model ' + str(hmm_group) + ' Search Hits</title>']
                html_report_lines += ['</head>']
                html_report_lines += ['<body bgcolor="white">']
                #if many_type_name == 'GenomeSet' or many_type_name == 'AnnotatedMetagenomeAssembly':
                html_report_lines += ['<a href="' + html_profile_file + '"><font color="' +
                                      header_tab_color + '" size=' + header_tab_fontsize + '>TABULAR PROFILE</font></a> | ']
                for this_hmm_group_i, this_hmm_group in enumerate(hmm_groups_used):
                    #disp_hmm_group = this_hmm_group[0].upper() + this_hmm_group[1:]
                    disp_hmm_group = model_group_config['fam_groups_disp'][this_hmm_group]
                    if this_hmm_group == hmm_group:
                        html_report_lines += [' <font color="' + header_tab_color + '" size=' +
                                              header_tab_fontsize + '><b>' + disp_hmm_group + ' HITS</b></font> ']
                    else:
                        this_html_search_file = search_tool_name + '_Search-' + \
                            str(this_hmm_group_i) + '-' + str(this_hmm_group) + '.html'
                        html_report_lines += [' <a href="' + this_html_search_file + '"><font color="' + header_tab_color +
                                              '" size=' + header_tab_fontsize + '>' + str(disp_hmm_group) + ' HITS</font></a> ']
                    if this_hmm_group_i < len(hmm_groups_used) - 1:
                        html_report_lines += [' | ']

                html_report_lines += ['<p>']
                html_report_lines += ['<table cellpadding=' + cellpadding +
                                      ' cellspacing = ' + cellspacing + ' border=' + border + '>']
                html_report_lines += ['<tr bgcolor="' + head_color + '">']
                html_report_lines += ['<td style="border-right:solid 2px ' + border_head_color + '; border-bottom:solid 2px ' + border_head_color +
                                      '"><font color="' + text_color + '" size=' + text_fontsize + '>' + 'ALIGNMENT COVERAGE (HIT SEQ)' + '</font></td>']
                html_report_lines += ['<td style="border-right:solid 2px ' + border_head_color + '; border-bottom:solid 2px ' +
                                      border_head_color + '"><font color="' + text_color + '" size=' + text_fontsize + '>' + 'GENE ID' + '</font></td>']
                html_report_lines += ['<td style="border-right:solid 2px ' + border_head_color + '; border-bottom:solid 2px ' +
                                      border_head_color + '"><font color="' + text_color + '" size=' + text_fontsize + '>' + 'FUNCTION' + '</font></td>']
                html_report_lines += ['<td style="border-right:solid 2px ' + border_head_color + '; border-bottom:solid 2px ' +
                                      border_head_color + '"><font color="' + text_color + '" size=' + text_fontsize + '>' + 'GENOME' + '</font></td>']
                #            html_report_lines += ['<td align=center style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'IDENT'+'%</font></td>']
                html_report_lines += ['<td align=center  style="border-right:solid 2px ' + border_head_color + '; border-bottom:solid 2px ' +
                                      border_head_color + '"><font color="' + text_color + '" size=' + text_fontsize + '>' + 'ALN_LEN' + '</font></td>']
                html_report_lines += ['<td align=center  style="border-right:solid 2px ' + border_head_color + '; border-bottom:solid 2px ' +
                                      border_head_color + '"><font color="' + text_color + '" size=' + text_fontsize + '>' + 'i-EVALUE' + '</font></td>']
                html_report_lines += ['<td align=center  style="border-right:solid 2px ' + border_head_color + '; border-bottom:solid 2px ' +
                                      border_head_color + '"><font color="' + text_color + '" size=' + text_fontsize + '>' + 'BIT SCORE' + '</font></td>']
                html_report_lines += ['<td align=center  style="border-right:solid 2px ' + border_head_color + '; border-bottom:solid 2px ' +
                                      border_head_color + '"><font color="' + text_color + '" size=' + text_fontsize + '>' + '<nobr>H_BEG-H_END</nobr>' + '</font></td>']
                #            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'MIS MATCH'+'</font></td>']
                #            html_report_lines += ['<td align=center  style="border-right:solid 2px '+border_head_color+'; border-bottom:solid 2px '+border_head_color+'"><font color="'+text_color+'" size='+text_fontsize+'>'+'GAP OPEN'+'</font></td>']
                html_report_lines += ['</tr>']

                for hmm_i, hmm_id in enumerate(input_HMM_ids[hmm_group]):
                    html_report_lines += ['<tr><td colspan=table_col_width>' +
                                          '<a name="'+str(hmm_id)+'"></a>' +
                                          'Hits to <b>' +str(hmm_id) +
                                          '</b></td></tr>']
                    if hmm_id not in total_hit_cnts:
                        html_report_lines += ['<tr><td colspan=table_col_width><blockquote><i>Model '+hmm_id+' not requested</i></td></tr>']
                    elif total_hit_cnts[hmm_id] == 0 or html_report_chunks[hmm_i] == None or html_report_chunks[hmm_i] == '':
                        html_report_lines += ['<tr><td colspan=table_col_width><blockquote><i>no hits found</i></td></tr>']
                    else:
                        #html_report_lines.extend(html_report_chunks[hmm_i])
                        html_report_lines += [html_report_chunks[hmm_i]]
                    html_report_lines += ['<tr><td colspan=table_col_width>' + sp + '</td></tr>']

                html_report_lines += ['</table>']
                html_report_lines += ['</body>']
                html_report_lines += ['</html>']

                # write html to file
                html_path = html_search_path
                html_report_str = "\n".join(html_report_lines)
                with open(html_path, 'w', 0) as html_handle:
                    html_handle.write(html_report_str)


        #### Build Profile output report
        ##
        self.log(console, "BUILDING PROFILE REPORT ")  # DEBUG
        #if len(invalid_msgs) == 0 and many_type_name == 'GenomeSet':
        if len(invalid_msgs) == 0:

            # calculate table
            #
            cats = all_HMM_ids_order
            table_data = dict()
            table_genes = dict()
            INSANE_VALUE = 10000000000000000
            if params.get('count_category') and params['count_category'] == 'perc_all':
                overall_low_val = INSANE_VALUE
            elif params.get('low_val') and params['low_val'] != 'detect':
                overall_low_val = float(params['low_val'])
            else:
                overall_low_val = INSANE_VALUE
            overall_high_val = -INSANE_VALUE
            cat_seen = dict()
            for cat in cats:
                cat_seen[cat] = False

            # count raw
            for genome_ref in all_genome_refs:
                if genome_ref not in table_data:
                    table_data[genome_ref] = dict()
                    table_genes[genome_ref] = dict()
                for cat in cats:
                    table_data[genome_ref][cat] = 0
                    table_genes[genome_ref][cat] = []

                if genome_ref not in hit_cnt_by_genome_and_model:
                    continue

                for cat in cats:
                    if cat in hit_cnt_by_genome_and_model[genome_ref] and \
                       hit_cnt_by_genome_and_model[genome_ref][cat] != 0:
                        table_data[genome_ref][cat] = hit_cnt_by_genome_and_model[genome_ref][cat]
                        table_genes[genome_ref][cat] = hit_genes_by_genome_and_model[genome_ref][cat]
                        cat_seen[cat] = True

            # adjust to perc all CDS if not raw count
            if params.get('count_category') and params['count_category'] == 'perc_all':
                for genome_ref in all_genome_refs:
                    total_genes = genome_CDS_count_by_ref[genome_ref]
                    for cat in cats:
                        if table_data[genome_ref][cat] != 0:
                            table_data[genome_ref][cat] /= float(total_genes)
                            table_data[genome_ref][cat] *= 100.0
                        
            # determine high and low val
            for genome_ref in all_genome_refs:
                for cat in cats:
                    val = table_data[genome_ref][cat]
                    if val == 0:
                        continue
                    #self.log (console, "HIGH VAL SCAN CAT: '"+cat+"' VAL: '"+str(val)+"'")  # DEBUG
                    if val > overall_high_val:
                        overall_high_val = val
                    if val < overall_low_val:
                        overall_low_val = val
            if overall_high_val == -INSANE_VALUE:
                error_msg = "unable to find any counts"
                self.log(invalid_msgs, error_msg)
                message = "\n".join(invalid_msgs)
                return self._create_error_report (params['workspace_name'], message, provenance)

            
            # build html report
            sp = '&nbsp;'
            body_bgcolor = '#ffffff'
            text_color = "#606060"
            text_color_2 = "#606060"
            head_color_1 = "#eeeeee"
            head_color_2 = "#eeeeee"
            border_color = "#cccccc"
            border_cat_color = "#ffccff"
            #graph_color = "lightblue"
            #graph_width = 100
            #graph_char = "."
            graph_char = sp
            #color_list = ['0','1','2','3','4','5','6','7','8','9','a','b','c','d','e']
            #color_list = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd']
            color_list = [
                "#333333",
                "#222266",
                "#222299",
                "#2222bb",
                "#2222dd",
                "#2222ff",
                "#4444ff",
                "#6666ff",
                "#8888ff",
                "#aaaaff",
                "#ccccff"]
            max_color = len(color_list) - 1
            cat_disp_trunc_len = 40
            cell_width = '20'
            cell_height = '18'
            #corner_radius = str(int(0.2*int(cell_width)+0.5))
            corner_radius = '5'
            if len(all_genome_refs) > 20:
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
            #graph_padding = "5"
            graph_padding = "1"
            graph_spacing = "3"
            #border = "1"
            border = "0"
            #row_spacing = "-2"
            num_rows = len(all_genome_refs)
            if int(params.get('show_target_block_headers',1)) == 1 and \
                len(input_many_refs) > 1:
                num_rows += len(input_many_refs)
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

            html_report_lines += [
                "/* unvisited link */\na:link {\n  color: "+text_color+";\n  text-decoration: none;\n}\n/* visited link */\na:visited {\n  color: "+text_color+";\n  text-decoration: none;\n}\n/* mouse over link */\na:hover {\n  color: blue;\n  text-decoration: underline;\n}\n/* selected link */\na:active {\n  color: red;\n  text-decoration: underline;\n}\n"]
            html_report_lines += [
                ".horz-text {\ndisplay: inline-block;\nfont-family: Tahoma, Geneva, sans-serif;\ntext-decoration: none;\n}"]
            html_report_lines += [
                ".vertical-text {\ndisplay: inline-block;\nfont-family: Tahoma, Geneva, sans-serif;\ntext-decoration: none;\nwidth: 0.65em;\n}\n.vertical-text__inner {\ndisplay: inline-block;\nwhite-space: nowrap;\nline-height: 1.1;\ntransform: translate(0,100%) rotate(-90deg);\ntransform-origin: 0 0;\n}\n.vertical-text__inner:after {\ncontent: \"\";\ndisplay: block;\nmargin: 0.0em 0 100%;\n}"]
            html_report_lines += [
                ".vertical-text_title {\ndisplay: inline-block;\nwidth: 1.0em;\n}\n.vertical-text__inner_title {\ndisplay: inline-block;\nwhite-space: nowrap;\nline-height: 1.0;\ntransform: translate(0,100%) rotate(-90deg);\ntransform-origin: 0 0;\n}\n.vertical-text__inner_title:after {\ncontent: \"\";\ndisplay: block;\nmargin: 0.0em 0 100%;\n}"]

            # add colors as style for DIV
            html_report_lines += [".heatmap_cell-"+'BLANK'+" {\nwidth: "+str(cell_width)+"px;\nheight: "+str(cell_height)+"px;\nborder-radius: "+str(corner_radius)+"px;\nbackground-color: "+body_bgcolor+";\ntext-align: center;\n}"]
            
            for color_i,color_val in enumerate(color_list):
                html_report_lines += [".heatmap_cell-"+str(color_i)+" {\nwidth: "+str(cell_width)+"px;\nheight: "+str(cell_height)+"px;\nborder-radius: "+str(corner_radius)+"px;\nbackground-color: "+str(color_val)+";\ntext-align: center;\n}"]
            
            html_report_lines += ['</style>']
            html_report_lines += ['</head>']
            html_report_lines += ['<body bgcolor="'+body_bgcolor+'">']
            html_report_lines += ['<font color="' + header_tab_color + '" size=' +
                                  header_tab_fontsize + '><b>TABULAR PROFILE</b></font> | ']

            for this_hmm_group_i, this_hmm_group in enumerate(hmm_groups_used):
                #disp_hmm_group = this_hmm_group[0].upper() + this_hmm_group[1:]
                disp_hmm_group = model_group_config['fam_groups_disp'][this_hmm_group]
                this_html_search_file = search_tool_name + '_Search-' + \
                    str(this_hmm_group_i) + '-' + str(this_hmm_group) + '.html'
                html_report_lines += [' <a href="' + this_html_search_file + '"><font color="' + header_tab_color +
                                      '" size=' + header_tab_fontsize + '>' + str(disp_hmm_group) + ' HITS</font></a> ']
                if this_hmm_group_i < len(hmm_groups_used) - 1:
                    html_report_lines += [' | ']
            html_report_lines += ['<p>']

            # genomes as rows
            if 'vertical' in params and int(params['vertical']) == 1:
                # table header
                html_report_lines += ['<table cellpadding=' + graph_padding +
                                      ' cellspacing=' + graph_spacing + ' border=' + border + '>']
                corner_rowspan = "2"
                label = ''
                html_report_lines += ['<tr>']
                html_report_lines += ['<td valign=bottom align=right rowspan=' + corner_rowspan +
                                      '><div class="vertical-text_title"><div class="vertical-text__inner_title"><font color="' + text_color + '">' + sp + label + '</font></div></div></td>']

                # fam groups
                for fam_group in fam_groups:
                    fam_seen_width = 0
                    for cat in all_HMM_ids[fam_group]:
                        if not cat_seen.get(cat) and not show_blanks:
                            continue
                        fam_seen_width += 1
                    if fam_seen_width == 0:
                        continue
                    fam_group_disp = model_group_config['fam_groups_disp'][fam_group]
                    cell_title = fam_group_disp
                    html_report_lines += ['<td colspan='+str(fam_seen_width)+' style="border-right:solid 2px ' + border_cat_color + '; border-bottom: solid 2px ' +
                                          border_cat_color + '" bgcolor="' + head_color_2 + '" title="' + cell_title + '" valign=bottom align=center>']
                    html_report_lines += ['<div class="horz_text">']
                    html_report_lines += ['<font color="'+text_color+'">']
                    html_report_lines += [fam_group_disp]
                    html_report_lines += ['</font.']
                    html_report_lines += ['</div>']
                    html_report_lines += ['</td>']
                html_report_lines += ['</tr>']
                        
                # column headers
                html_report_lines += ['<tr>']
                for cat_i, cat in enumerate(cats):
                    if not cat_seen.get(cat) and not show_blanks:
                        continue
                    cat_disp = cat
                    #cell_title = input_HMM_descs[cat]
                    cell_title = cat
                    key_link = cat
                    
                    if len(cat_disp) > cat_disp_trunc_len + 1:
                        cat_disp = cat_disp[0:cat_disp_trunc_len] + '*'

                    html_report_lines += ['<td style="border-right:solid 2px ' + border_cat_color + '; border-bottom: solid 2px ' +
                                          border_cat_color + '" bgcolor="' + head_color_2 + '" title="' + cell_title + '" valign=bottom align=center>']
                    html_report_lines += ['<div class="vertical-text"><div class="vertical-text__inner">']
                    html_report_lines += ['<font color="' + text_color_2 + '" size=' + graph_cat_fontsize + '><b>']
                    
                    #for c_i,c in enumerate(cat_disp):
                    #    if c_i < len(cat_disp)-1:
                    #        html_report_lines += [c+'<br>']
                    #    else:
                    #        html_report_lines += [c]
                    html_report_lines += ['<nobr>'+sp]
                    html_report_lines += ['<a href="#'+key_link+'" target="_key_tab">']
                    html_report_lines += [cat_disp]
                    html_report_lines += ['</a></nobr>']
                    html_report_lines += ['</b></font>']
                    html_report_lines += ['</div></div>']
                    html_report_lines += ['</td>']
                html_report_lines += ['</tr>']


                # rest of rows
                num_cols = 0;
                for cat in cats:
                    if not cat_seen.get(cat) and not show_blanks:
                        continue
                    num_cols += 1

                for input_many_ref in input_many_refs:
                    many_type_name = many_type_names[input_many_ref]

                    if int(params.get('show_target_block_headers',1)) == 1 and \
                        len(input_many_refs) > 1:
                        input_obj_disp_name = input_many_names[input_many_ref]
                        html_report_lines += ['<tr>']
                        html_report_lines += ['<td align=right><div class="horz-text"><font color="' + text_color + '" size=' +
                                          graph_gen_fontsize + '><b><nobr>' + input_obj_disp_name + sp + '</nobr></b></font></div></td>']
                        html_report_lines += ['<td colspan='+str(num_cols)+' bgcolor="#dddddd"></td>']
                        html_report_lines += ['</tr>']
                    

                    # add genome and ama rows, order GenomeSet elements by disp name
                    genome_ref_order = []
                    if many_type_name != 'GenomeSet' and many_type_name != 'FeatureSet' and many_type_name != 'Tree':
                        genome_ref_order = genome_refs[input_many_ref]
                        
                    # Use tree order, first ETE3 ladderize tree
                    elif many_type_name == 'Tree':
                        input_many_tree = input_many_objs[input_many_ref]
                        genome_id_by_ref = dict()
                        genome_ref_by_id = dict()    
                        for genome_id in input_many_tree['default_node_labels'].keys():
                            genome_ref = input_many_tree['ws_refs'][genome_id]['g'][0]
                            genome_id_by_ref[genome_ref] = genome_id
                            genome_ref_by_id[genome_id] = genome_ref
                        species_tree = ete3.Tree(input_many_tree['tree'])
                        species_tree.ladderize()
                        for genome_id in species_tree.get_leaf_names():
                            genome_ref_order.append(genome_ref_by_id[genome_id])

                    # Sort disp names alphabetically
                    else:
                        genome_disp_name_to_ref = dict()
                        for genome_ref in genome_refs[input_many_ref]:
                            genome_disp_name = self._get_genome_disp_name (params,
                                                                           genome_ref,
                                                                           many_type_name,
                                                                           ama_ref_to_obj_name,
                                                                           genome_ref_to_obj_name,
                                                                           genome_ref_to_sci_name)
                            if genome_disp_name not in genome_disp_name_to_ref:
                                genome_disp_name_to_ref[genome_disp_name] = []
                            genome_disp_name_to_ref[genome_disp_name].append(genome_ref)

                        for genome_disp_name in sorted(genome_disp_name_to_ref.keys()):
                            for genome_ref in genome_disp_name_to_ref[genome_disp_name]:
                                genome_ref_order.append(genome_ref)

                    for genome_ref in genome_ref_order:
                        genome_disp_name = self._get_genome_disp_name (params,
                                                                       genome_ref,
                                                                       many_type_name,
                                                                       ama_ref_to_obj_name,
                                                                       genome_ref_to_obj_name,
                                                                       genome_ref_to_sci_name)

                        # build html report table line
                        html_report_lines += ['<tr>']
                        html_report_lines += ['<td align=right><div class="horz-text"><font color="' + text_color + '" size=' +
                                              graph_gen_fontsize + '><b><nobr>' + genome_disp_name + sp + '</nobr></b></font></div></td>']
                        for cat in cats:
                            if not cat_seen.get(cat) and not show_blanks:
                                continue
                            val = table_data[genome_ref][cat]
                            if not cat_seen.get(cat) or val == 0:
                                #html_report_lines += ['<td bgcolor=white></td>']
                                html_report_lines += ['<td align=center valign=middle bgcolor="'+body_bgcolor+'">']
                                html_report_lines += ['<div class="heatmap_cell-BLANK"></div>']
                                html_report_lines += ['</td>']
                                continue
                            elif overall_high_val == overall_low_val:
                                cell_color_i = 0
                            else:
                                cell_color_i = max_color - \
                                               int(round(max_color * (val - overall_low_val) / float(overall_high_val - overall_low_val)))

                            cell_val = str(table_data[genome_ref][cat])  # the key line

                            if 'heatmap' in params and params['heatmap'] == '1':
                                s = 's'
                                if cell_val == '1': s = ''
                                cell_title = cell_val+' hit'+s+"\n"+"\n".join(sorted(table_genes[genome_ref][cat]))
	                        hmm_group = hmm_id_to_hmm_group_name[cat]
                                hmm_group_i = hmm_id_to_hmm_group_index[cat]
                                group_html_search_file = search_tool_name + '_Search-' + \
                                                         str(hmm_group_i) + '-' + str(hmm_group) + '.html'
                                html_report_lines += ['<td title="'+cell_title+'" align=center valign=middle bgcolor="'+body_bgcolor+'">' +
                                                      '<a href="'+group_html_search_file+'#'+cat+'">' +
                                                      '<div class="heatmap_cell-'+str(cell_color_i)+'"></div>' +
                                                      '</a></td>']
                            else:
                                html_report_lines += ['<td align=center valign=middle style="' + cell_width + 'px; border-right:solid 2px ' + border_color +
                                                      '; border-bottom:solid 2px ' + border_color + '"><font color="' + text_color + '" size=' + cell_fontsize + '>' + cell_val + '</font></td>']
                        
                        html_report_lines += ['</tr>']
                html_report_lines += ['</table>']

            # genomes as columns
            else:
                raise ValueError("Do not yet support Genomes as columns")

            # key table
            html_report_lines += ['<p>']
            html_report_lines += ['<table cellpadding=3 cellspacing=2 border=' + border + '>']
            html_report_lines += ['<tr><td valign=middle align=left colspan=2 style="border-bottom:solid 4px ' +
                                  border_color + '"><font color="' + text_color + '"><b>KEY</b></font></td></tr>']

            # fam groups
            for fam_group in fam_groups:
                fam_group_disp = model_group_config['fam_groups_disp'][fam_group]
                fam_group_col_width = '2'
                fam_group_cell_color = 'white'
                fam_group_text_color = text_color
                html_report_lines += ['<tr>']
                html_report_lines += ['<td valign=middle align=left bgcolor="' + fam_group_cell_color + '">']
                html_report_lines += ['<br>']
                html_report_lines += ['<font color="'+fam_group_text_color+'"><b>']
                html_report_lines += [fam_group_disp]
                html_report_lines += ['</b></font>']
                html_report_lines += ['</td>']
                html_report_lines += ['</tr>']

                for cat in all_HMM_ids[fam_group]:
                    cell_color = 'white'
                    #if not cat_seen.get(cat) and not show_blanks:
                    if not cat_seen.get(cat):
                        cell_color = "#eeeeee"
                    cat_desc = input_HMM_descs[cat]
                    key_link = cat
                    cat_disp = cat
                    if len(cat_disp) > cat_disp_trunc_len + 1:
                        cat_disp = cat_disp[0:cat_disp_trunc_len] + '*'

                    html_report_lines += ['<tr>']
                    html_report_lines += ['<td valign=middle align=left bgcolor="' + cell_color + '" style="border-right:solid 4px ' + border_color +
                                          '">' + '<div class="horz-text">' + '<font color="' + text_color + '" size=' + graph_cat_fontsize + '>' + cat_disp + '</font>' + '</div>' + '</td>']
                    html_report_lines += ['<td valign=middle align=left bgcolor="' + cell_color + '">' +
                                          '<div class="horz-text">' +
                                          '<a name="'+key_link+'">' +
                                          '<font color="' + text_color + '" size=' + graph_cat_fontsize + '>' + cat_desc + '</font>' + '</div>' + '</td>']
                    html_report_lines += ['</tr>']
                
            # close
            html_report_lines += ['</table>']
            html_report_lines += ['</body>']
            html_report_lines += ['</html>']

            # write html to file and upload
            html_path = html_profile_path
            html_report_str = "\n".join(html_report_lines)
            with open(html_path, 'w', 0) as html_handle:
                html_handle.write(html_report_str)

        #### Upload HTML reports
        ##
        self.log(console, "UPLOADING HTML REPORT(s)")  # DEBUG
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
                raise ValueError('Logging exception loading HTML file to shock')

        #### Upload output files
        ##
        self.log(console, "UPLOADING OUTPUT FILES")  # DEBUG
        if len(invalid_msgs) == 0:

            output_hit_TAB_dir = os.path.join(self.output_dir, 'HMMER_output_TAB')
            #output_hit_MSA_dir = os.path.join(self.output_dir, 'HMMER_output_MSA')
            if not os.path.exists(output_hit_TAB_dir):
                os.makedirs(output_hit_TAB_dir)
            #if not os.path.exists(output_hit_MSA_dir):
            #    os.makedirs(output_hit_MSA_dir)

            for hmm_i, hmm_id in enumerate(all_HMM_ids_order):
                if hmm_id not in total_hit_cnts:
                    self.log(console, 'MODEL NOT REQUESTED.  SKIPPING UPLOAD FOR HMM ' + hmm_id)
                    continue
                elif total_hit_cnts[hmm_id] == 0:
                    self.log(console, 'SKIPPING UPLOAD OF EMPTY HMMER OUTPUT FOR HMM ' + hmm_id)
                    continue
                else:
                    self.log(console, 'PREPPING UPLOAD OF HMMER OUTPUT FOR HMM ' + hmm_id)
                for input_many_ref in input_many_refs:
                    input_many_name = input_many_names[input_many_ref]
                    new_hit_TAB_file_path = os.path.join(output_hit_TAB_dir, hmm_id +'-'+ input_many_name + '.hitout.txt')
                    #new_hit_MSA_file_path = os.path.join(output_hit_MSA_dir, hmm_id + '.msaout.txt')

                    shutil.copy(output_hit_TAB_file_paths[hmm_id][input_many_ref], new_hit_TAB_file_path)
                    #shutil.copy(output_hit_MSA_file_paths[hmm_id], new_hit_MSA_file_path)

            # Upload output dirs
            TAB_upload_ret = None
            MSA_upload_ret = None
            self.log(console, 'UPLOADING OF HMMER OUTPUT')
            try:
                TAB_upload_ret = dfu.file_to_shock({'file_path': output_hit_TAB_dir,
                                                    'make_handle': 0,
                                                    'pack': 'zip'})
            except:
                raise ValueError('Logging exception loading TAB output to shock')
            '''
            try:
                MSA_upload_ret = dfu.file_to_shock({'file_path': output_hit_MSA_dir,
                                                    'make_handle': 0,
                                                    'pack': 'zip'})
            except:
                raise ValueError('Logging exception loading MSA output to shock')
            '''
            
        #### Create report object
        ##
        self.log(console, "CREATING REPORT OBJECT")  # DEBUG
        if len(invalid_msgs) == 0:

            reportName = 'hmmer_report_' + str(uuid.uuid4())
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
                                        'name': html_profile_file,
                                        'label': search_tool_name + ' HTML Report'}
                                        #'description': search_tool_name + ' HTML Report'}
                                       ]

            if TAB_upload_ret != None:
                reportObj['file_links'] += [{'shock_id': TAB_upload_ret['shock_id'],
                                             'name': search_tool_name + '_Search.TAB.zip',
                                             'label': search_tool_name + '-' + ' hits TABLE'}]
            if MSA_upload_ret != None:
                reportObj['file_links'] += [{'shock_id': MSA_upload_ret['shock_id'],

                                             'name': search_tool_name + '_Search.MSA.zip',
                                             'label': search_tool_name + ' hits MSA'}
                                            ]

            for hmm_group in all_HMM_groups_order:
                if hit_accept_something[hmm_group]:
                    if 'coalesce_output' in params and int(params['coalesce_output']) == 1:
                        if hmm_group in objects_created_refs_coalesce:
                            reportObj['objects_created'].append(
                                {'ref': objects_created_refs_coalesce[hmm_group], 'description': 'Coalesced' + ' ' + hmm_group + ' ' + search_tool_name + ' hits'})
                    else:
                        for hmm_i, hmm_id in enumerate(all_HMM_ids[hmm_group]):
                            if hmm_id in objects_created_refs_by_hmm_id:
                                reportObj['objects_created'].append(
                                    {'ref': objects_created_refs_by_hmm_id[hmm_id], 'description': hmm_id + ' ' + search_tool_name + ' hits'})

            # save report object
            #
            SERVICE_VER = 'release'
            reportClient = KBaseReport(self.callbackURL, token=ctx['token'], service_ver=SERVICE_VER)
            #report_info = report.create({'report':reportObj, 'workspace_name':params['workspace_name']})
            report_info = reportClient.create_extended_report(reportObj)

        #### data validation error
        ##
        if len(invalid_msgs) > 0:
            message = "\n".join(invalid_msgs)
            return self._create_error_report (params['workspace_name'], message, provenance)

        
        #### Return Report
        ##
        self.log(console, "BUILDING RETURN OBJECT")
#        returnVal = { 'output_report_name': reportName,
#                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
#                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
#                      }
        returnVal = {'report_name': report_info['name'],
                     'report_ref': report_info['ref']
                     }
        self.log(console, search_tool_name + "_Search DONE")

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method HMMER_Model_Group_Search return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return returnVal
