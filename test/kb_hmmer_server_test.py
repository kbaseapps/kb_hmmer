# -*- coding: utf-8 -*-
import unittest
import os
import json
import time
import shutil
import requests
import uuid

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint

from biokbase.workspace.client import Workspace as workspaceService
from SetAPI.SetAPIServiceClient import SetAPI
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil

from kb_hmmer.kb_hmmerImpl import kb_hmmer
from kb_hmmer.kb_hmmerServer import MethodContext
from kb_hmmer.authclient import KBaseAuth as _KBaseAuth


class kb_hmmerTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_hmmer'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_hmmer',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = kb_hmmer(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

        suffix = int(time.time() * 1000)
        cls.wsName = "test_kb_Msuite_" + str(suffix)
        cls.ws_info = cls.wsClient.create_workspace({'workspace': cls.wsName})
        cls.setAPI = SetAPI(url=cls.cfg['service-wizard-url'], token=cls.ctx['token'])
        cls.gfu = GenomeFileUtil(os.environ['SDK_CALLBACK_URL'])

        cls.prepare_data() # prepare WS data


    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_hmmer_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx


    @classmethod
    def prepare_data(cls):
        test_directory_name = 'test_kb_Msuite'
        cls.test_directory_path = os.path.join(cls.scratch, test_directory_name)
        os.makedirs(cls.test_directory_path)

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple


        # Upload a few genomes
        cls.genome_refs = []
#        for i,genome_filename in enumerate(['GCF_000018425.1_ASM1842v1_genomic.gbff']):  # DEBUG
        for i,genome_filename in enumerate(['GCF_000018425.1_ASM1842v1_genomic.gbff', \
                                            'GCF_000022285.1_ASM2228v1_genomic.gbff', \
                                            'GCF_001439985.1_wTPRE_1.0_genomic.gbff']): 
                                            #'GCF_000287295.1_ASM28729v1_genomic.gbff', \
                                            #'GCF_000306885.1_ASM30688v1_genomic.gbff', \

            print ("prepare_data(): UPLOADING GENOME: "+genome_filename)  # DEBUG
            genome_file_path = os.path.join(cls.scratch, genome_filename)
            shutil.copy(os.path.join("data", "genomes", genome_filename), genome_file_path)
            cls.genome_refs.append(cls.gfu.genbank_to_genome({'file': {'path': genome_file_path},
                                                              'workspace_name': cls.ws_info[NAME_I],
                                                              'genome_name': genome_filename})['genome_ref'])


        # Create a GenomeSet
        cls.genomeSet_refs = []
        genomeSet_name = 'test_genomeset_1'
        genome_scinames = dict()
        for genome_i,genome_ref in enumerate(cls.genome_refs):
            genome_scinames[genome_ref] = 'Genus species str. '+str(genome_i)
        testGS = {
            'description': 'genomeSet for testing',
            'elements': dict()
        }
        for genome_ref in cls.genome_refs: 
            testGS['elements'][genome_scinames[genome_ref]] = { 'ref': genome_ref }
        print ("prepare_data(): UPLOADING GENOME SET: "+genomeSet_name)  # DEBUG
        obj_info = cls.wsClient.save_objects({'workspace': cls.ws_info[NAME_I],       
                                              'objects': [
                                                  {
                                                      'type':'KBaseSearch.GenomeSet',
                                                      'data':testGS,
                                                      'name':genomeSet_name,
                                                      'meta':{},
                                                      'provenance':[
                                                          {
                                                              'service':'kb_Msuite',
                                                              'method':'test_CheckM'
                                                          }
                                                      ]
                                                  }]
                                          })[0]
        cls.genomeSet_refs.append(str(obj_info[WSID_I]) +'/'+ str(obj_info[OBJID_I]) +'/'+ str(obj_info[VERSION_I]))


        # Upload MSAs
        cls.MSA_refs = []
        MSA_data_dir = os.path.join('data','MSA')
        MSA_ext = '.json'
        MSA_obj_type = 'KBaseTrees.MSA'
        for (dirpath, dirnames, filenames) in os.walk(MSA_data_dir):
            for f in sorted(filenames):   # filenames.sort() does not work?!?!?!?!
                if not f.endswith(MSA_ext):
                    continue
                this_MSA_json_path = os.path.join(MSA_data_dir, f)
                this_MSA_name = f.rstrip(MSA_ext)

                with open (this_MSA_json_path, 'r', 0) as MSA_json_fh:
                    this_MSA_obj = json.load(MSA_json_fh)

                print ("prepare_data(): UPLOADING MSA: "+this_MSA_name)  # DEBUG
                provenance = [{}]
                MSA_info = cls.wsClient.save_objects({
                    'workspace': cls.ws_info[NAME_I],
                    'objects': [
                        {
                            'type': MSA_obj_type,
                            'data': this_MSA_obj,
                            'name': this_MSA_name,
                            'meta': {},
                            'provenance': provenance
                        }
                    ]})[0]
                cls.MSA_refs.append(str(MSA_info[WSID_I])+'/'+str(MSA_info[OBJID_I])+'/'+str(MSA_info[VERSION_I]))



    #
    # NOTE: According to Python unittest naming rules test method names should start from 'test'.
    #


    ### Test 01: Single Model against Single Genome
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_01_kb_hmmer_HMMER_MSA_Search_Genome()")
    def test_01_kb_hmmer_HMMER_MSA_Search_Genome(self):
        test_name = 'test_01_kb_hmmer_HMMER_MSA_Search_Genome'
        header_msg = "RUNNING "+test_name+"()"
        header_delim = len(header_msg) * '='
        print ("\n"+header_delim+"\n"+header_msg+"\n"+header_delim+"\n")

        obj_basename = test_name+'.HMMER_MSA'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'

        # app run params
        parameters = { 'workspace_name': self.getWsName(),
                       'input_msa_ref': self.MSA_refs[0],      # Single MSA
                       'input_many_ref': self.genome_refs[0],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'e_value': ".001",
                       'bitscore': "50",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000"
                     }
        ret = self.getImpl().HMMER_MSA_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)
        pass


    ### Test 02: Single Model against GenomeSet
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_02_kb_hmmer_HMMER_MSA_Search_GenomeSet()")
    def test_02_kb_hmmer_HMMER_MSA_Search_GenomeSet(self):
        test_name = 'test_02_kb_hmmer_HMMER_MSA_Search_GenomeSet'
        header_msg = "RUNNING "+test_name+"()"
        header_delim = len(header_msg) * '='
        print ("\n"+header_delim+"\n"+header_msg+"\n"+header_delim+"\n")

        obj_basename = test_name+'.HMMER_MSA'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'

        # app run params        
        parameters = { 'workspace_name': self.getWsName(),
                       'input_msa_ref': self.MSA_refs[0],         # Single MSA
                       'input_many_ref': self.genomeSet_refs[0],  # GenomeSet
                       'output_filtered_name': obj_out_name,
                       'e_value': ".001",
                       'bitscore': "50",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000"
                     }
        ret = self.getImpl().HMMER_MSA_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)
        pass


    ### Test 03: All Models in workspace against Single Genome
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_03_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome()")
    def test_03_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome(self):
        test_name = 'test_03_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome'
        header_msg = "RUNNING "+test_name+"()"
        header_delim = len(header_msg) * '='
        print ("\n"+header_delim+"\n"+header_msg+"\n"+header_delim+"\n")

        obj_basename = test_name+'.HMMER_MSA'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'

        # app run params
        parameters = { 'workspace_name': self.getWsName(),
                       'use_all_local_MSAs': "0",
                       'input_msa_refs': [self.MSA_refs[0], self.MSA_refs[1], self.MSA_refs[2]],  # Specific MSAs
                       'input_many_ref': self.genome_refs[0],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'coalesce_output': 0,
                       'e_value': ".001",
                       'bitscore': "50",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_Local_MSA_Group_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_objs_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})
        for created_obj_info in created_objs_info:
            #self.assertEqual(created_obj_info[NAME_I], obj_out_name)  # MSA name is prepended
            self.assertEqual(created_obj_info[TYPE_I].split('-')[0], obj_out_type)
        pass


    ### Test 04: All Models in workspace against GenomeSet, DON'T coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_04_kb_hmmer_HMMER_Local_MSA_Group_Search_GenomeSet_NOcoalesce()")
    def test_04_kb_hmmer_HMMER_Local_MSA_Group_Search_GenomeSet_NOcoalesce(self):
        test_name = 'test_04_kb_hmmer_HMMER_Local_MSA_Group_Search_GenomeSet_NOcoalesce'
        header_msg = "RUNNING "+test_name+"()"
        header_delim = len(header_msg) * '='
        print ("\n"+header_delim+"\n"+header_msg+"\n"+header_delim+"\n")

        obj_basename = test_name+'.HMMER_MSA'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'

        # app run params        
        parameters = { 'workspace_name': self.getWsName(),
                       'use_all_local_MSAs': "1",
                       #'input_msa_refs': [self.MSA_refs[0]],         # Single MSA
                       'input_many_ref': self.genomeSet_refs[0],  # GenomeSet
                       'output_filtered_name': obj_out_name,
                       'coalesce_output': 0,  # KEY
                       'e_value': ".001",
                       'bitscore': "50",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_Local_MSA_Group_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created objs
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_objs_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})
        for created_obj_info in created_objs_info:
            #self.assertEqual(created_obj_info[NAME_I], obj_out_name)  # MSA name is prepended
            self.assertEqual(created_obj_info[TYPE_I].split('-')[0], obj_out_type)
        pass


    ### Test 05: All Models in workspace against GenomeSet, DO coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_05_kb_hmmer_HMMER_Local_MSA_Group_Search_GenomeSet_coalesce()")
    def test_05_kb_hmmer_HMMER_Local_MSA_Group_Search_GenomeSet_coalesce(self):
        test_name = 'test_05_kb_hmmer_HMMER_Local_MSA_Group_Search_GenomeSet_coalesce'
        header_msg = "RUNNING "+test_name+"()"
        header_delim = len(header_msg) * '='
        print ("\n"+header_delim+"\n"+header_msg+"\n"+header_delim+"\n")

        obj_basename = test_name+'.HMMER_MSA'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'

        # app run params        
        parameters = { 'workspace_name': self.getWsName(),
                       'use_all_local_MSAs': "1",
                       #'input_msa_ref': [self.MSA_refs[0]],         # Single MSA
                       'input_many_ref': self.genomeSet_refs[0],  # GenomeSet
                       'output_filtered_name': obj_out_name,
                       'coalesce_output': 1,  # KEY
                       'e_value': ".001",
                       'bitscore': "50",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_Local_MSA_Group_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created objs
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_objs_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})
        for created_obj_info in created_objs_info:
            #self.assertEqual(created_obj_info[NAME_I], obj_out_name)  # MSA name is prepended
            self.assertEqual(created_obj_info[TYPE_I].split('-')[0], obj_out_type)
        pass


    ### Test 06: Single Model against Single Genome, threshold above all hits
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_06_kb_hmmer_HMMER_MSA_Search_Genome_removeALL()")
    def test_06_kb_hmmer_HMMER_MSA_Search_Genome_removeALL(self):
        test_name = 'test_06_kb_hmmer_HMMER_MSA_Search_Genome_removeALL'
        header_msg = "RUNNING "+test_name+"()"
        header_delim = len(header_msg) * '='
        print ("\n"+header_delim+"\n"+header_msg+"\n"+header_delim+"\n")

        obj_basename = test_name+'.HMMER_MSA'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'

        # app run params
        parameters = { 'workspace_name': self.getWsName(),
                       'input_msa_ref': self.MSA_refs[0],      # Single MSA
                       'input_many_ref': self.genome_refs[0],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'e_value': ".001",
                       'bitscore': "5000000",
                       'overlap_fraction': "100.0",
                       'maxaccepts': "1000"
                     }
        ret = self.getImpl().HMMER_MSA_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_obj_0_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})[0]
        self.assertEqual(created_obj_0_info[NAME_I], obj_out_name)
        self.assertEqual(created_obj_0_info[TYPE_I].split('-')[0], obj_out_type)
        pass


    ### Test 07: All Models in workspace against Single Genome, threshold above all hits
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_07_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome_removeALL()")
    def test_07_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome_removeALL(self):
        test_name = 'test_07_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome_removeALL'
        header_msg = "RUNNING "+test_name+"()"
        header_delim = len(header_msg) * '='
        print ("\n"+header_delim+"\n"+header_msg+"\n"+header_delim+"\n")

        obj_basename = test_name+'.HMMER_MSA'
        obj_out_name = obj_basename+".test_output.FS"
        obj_out_type = "KBaseCollections.FeatureSet"

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        #reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        #genome_ref_1 = 'ReferenceDataManager/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'

        # app run params
        parameters = { 'workspace_name': self.getWsName(),
                       'use_all_local_MSAs': "0",
                       'input_msa_refs': [self.MSA_refs[0], self.MSA_refs[1], self.MSA_refs[2]],  # Specific MSAs
                       'input_many_ref': self.genome_refs[0],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'coalesce_output': 0,
                       'e_value': ".001",
                       'bitscore': "500000000",
                       'overlap_fraction': "100.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_Local_MSA_Group_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])

        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertIsNotNone(report_obj['objects_created'][0]['ref'])

        created_objs_info = self.getWsClient().get_object_info_new({'objects':[{'ref':report_obj['objects_created'][0]['ref']}]})
        for created_obj_info in created_objs_info:
            #self.assertEqual(created_obj_info[NAME_I], obj_out_name)  # MSA name is prepended
            self.assertEqual(created_obj_info[TYPE_I].split('-')[0], obj_out_type)
        pass
