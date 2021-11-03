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

from installed_clients.WorkspaceClient import Workspace
from installed_clients.SetAPIServiceClient import SetAPI
from installed_clients.GenomeFileUtilClient import GenomeFileUtil

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
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_hmmer(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

        suffix = int(time.time() * 1000)
        cls.wsName = "test_kb_hmmer_" + str(suffix)
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
        test_directory_name = 'test_kb_hmmer'
        cls.test_directory_path = os.path.join(cls.scratch, test_directory_name)
        os.makedirs(cls.test_directory_path)

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        # Upload a few genomes
        cls.genome_refs = []
        for i,genome_filename in enumerate(['GCF_000018425.1_ASM1842v1_genomic.gbff', \
                                            'GCF_000022285.1_ASM2228v1_genomic.gbff', \
                                            'GCF_001439985.1_wTPRE_1.0_genomic.gbff', \
                                            'GCF_000015865.1_ASM1586v1_genomic.gbff', \
                                            'GCF_000164865.1_ASM16486v1_genomic.gbff']):
                                            #'GCF_000287295.1_ASM28729v1_genomic.gbff', \
                                            #'GCF_000306885.1_ASM30688v1_genomic.gbff', \
        #for i,genome_filename in enumerate(['GCF_000018425.1_ASM1842v1_genomic.gbff']):  # DEBUG
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
                                                              'service':'kb_hmmer',
                                                              'method':'test_hmmer'
                                                          }
                                                      ]
                                                  }]
                                          })[0]
        cls.genomeSet_refs.append(str(obj_info[WSID_I]) +'/'+ str(obj_info[OBJID_I]) +'/'+ str(obj_info[VERSION_I]))


        # Upload some test Annotated Metagenome Assemblies
        cls.ama_refs = []
        ama_name = "ama_test.AMA"
        ama_feature_cnt = 888
        ama_contigs_file_src = "data/AnnotatedMetagenomeAssembly/ama_contigs.fasta"
        ama_genes_file_src   = "data/AnnotatedMetagenomeAssembly/ama_genes.gff"
        #ama_contigs_file_src = "data/AnnotatedMetagenomeAssembly/Unbinned_plus_LQ_bins.fasta"
        #ama_genes_file_src   = "data/AnnotatedMetagenomeAssembly/Unbinned_plus_LQ_bins.gff"
        ama_contigs_file_upload = os.path.join (cls.scratch, os.path.basename(ama_contigs_file_src))
        ama_genes_file_upload = os.path.join (cls.scratch, os.path.basename(ama_genes_file_src))
        shutil.copy (ama_contigs_file_src, ama_contigs_file_upload)
        shutil.copy (ama_genes_file_src, ama_genes_file_upload)
        ama_upload_params = {
            "workspace_name": cls.wsName,
            "fasta_file": {"path": ama_contigs_file_upload},
            "gff_file": {"path": ama_genes_file_upload},
            "source": "GFF",
            "scientific_name": "TEST AMA",
            "generate_missing_genes": "True"
        }
        try:
            SERVICE_VER = 'dev'
            GFU = GenomeFileUtil(os.environ['SDK_CALLBACK_URL'],
                                 token=cls.ctx['token'],
                                 service_ver=SERVICE_VER
                             )
        except:
            raise ValueError("unable to instantiate GFU client")
        for i in range(2):
            ama_upload_params['genome_name'] = ama_name+'-'+str(i)
            print ("UPLOADING AMA: "+ama_name+" to WORKSPACE "+cls.wsName+" ...")
            try:
                ama_upload_result = GFU.fasta_gff_to_metagenome (ama_upload_params)
            except:
                raise ValueError("unable to upload test AMA data object")
            pprint (ama_upload_result)
            ama_ref = ama_upload_result['metagenome_ref']
            cls.ama_refs.append(ama_ref)
        

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
                       'input_msa_ref': self.MSA_refs[1],      # Single MSA
                       'input_many_ref': self.genome_refs[0],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name',
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
                       'input_msa_ref': self.MSA_refs[1],         # Single MSA
                       'input_many_ref': self.genomeSet_refs[0],  # GenomeSet
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'sci_name',
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
    
    
    ### Test 03: Single nucleotide Model (to test failure)
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_03_kb_hmmer_HMMER_MSA_Search_Genome_nuc_MSA()")
    def test_03_kb_hmmer_HMMER_MSA_Search_Genome_nuc_MSA(self):
        test_name = 'test_03_kb_hmmer_HMMER_MSA_Search_Genome_nuc_MSA'
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
                       'genome_disp_name_config': 'obj_name_sci_name',
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
        self.assertEqual(report_obj['text_message'][0:7],"FAILURE")
        pass
    
    
    ### Test 03_02: Single Model against Single AnnotatedMetagenomeAssembly
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_03_02_kb_hmmer_HMMER_MSA_Search_AMA()")
    def test_003_02_kb_hmmer_HMMER_MSA_Search_AMA(self):
        test_name = 'test_03_02_kb_hmmer_HMMER_MSA_Search_AMA'
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
                       'input_msa_ref': self.MSA_refs[1],      # Single MSA
                       'input_many_ref': self.ama_refs[0],  # Single AnnotatedMetagenomeAssembly
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name',
                       #'e_value': ".001",
                       'e_value': ".1",
                       #'bitscore': "50",
                       'bitscore': "1",
                       #'overlap_fraction': "50.0",
                       'overlap_fraction': "5.0",
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
    
    
    ### Test 04: Specific Models in workspace against Single Genome
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_04_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome()")
    def test_04_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome(self):
        test_name = 'test_04_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome'
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
                       'input_msa_refs': [self.MSA_refs[1], self.MSA_refs[2], self.MSA_refs[3]],  # Specific MSAs
                       'input_many_ref': self.genome_refs[0],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,
                       'e_value': ".001",
                       'bitscore': "50",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "1",
                       'vertical': "1",
                       'show_blanks': "1"
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
    
    
    ### Test 04_02: Specific Models in workspace against Single Annotated Metagenome Assembly
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_04_02_kb_hmmer_HMMER_Local_MSA_Group_Search_AMA()")
    def test_04_02_kb_hmmer_HMMER_Local_MSA_Group_Search_AMA(self):
        test_name = 'test_04_02_kb_hmmer_HMMER_Local_MSA_Group_Search_AMA'
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
                       'input_msa_refs': [self.MSA_refs[1], self.MSA_refs[2], self.MSA_refs[3]],  # Specific MSAs
                       'input_many_ref': self.ama_refs[0],  # Single AnnotatedMetagenomeAssembly
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,
                       #'e_value': ".001",
                       'e_value': ".1",
                       #'bitscore': "50",
                       'bitscore': "10",
                       #'overlap_fraction': "50.0",
                       'overlap_fraction': "5.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "0.1",
                       'vertical': "1",
                       'show_blanks': "1"
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


    ### Test 04_03: All Models in workspace against Single Annotated Metagenome Assembly, DON'T Coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_04_04_kb_hmmer_HMMER_Local_MSA_Group_Search_AMA_ALL_DONT_coalesce()")
    def test_04_03_kb_hmmer_HMMER_Local_MSA_Group_Search_AMA_ALL_DONT_coalesce(self):
        test_name = 'test_04_03_kb_hmmer_HMMER_Local_MSA_Group_Search_AMA_ALL_DONT_coalesce'
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
                       'input_msa_refs': [],
                       'input_many_ref': self.ama_refs[0],  # Single AnnotatedMetagenomeAssembly
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,
                       #'e_value': ".001",
                       'e_value': ".1",
                       #'bitscore': "50",
                       'bitscore': "10",
                       #'overlap_fraction': "50.0",
                       'overlap_fraction': "5.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "detect",
                       'vertical': "1",
                       'show_blanks': "1"
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

    
    ### Test 05: All Models in workspace against GenomeSet, DON'T coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_05_kb_hmmer_HMMER_Local_MSA_Group_Search_GenomeSet_NOcoalesce()")
    def test_05_kb_hmmer_HMMER_Local_MSA_Group_Search_GenomeSet_NOcoalesce(self):
        test_name = 'test_05_kb_hmmer_HMMER_Local_MSA_Group_Search_GenomeSet_NOcoalesce'
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
                       'input_many_ref': self.genomeSet_refs[0],  # GenomeSet
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_ver_sci_name',
                       'coalesce_output': 0,  # KEY
                       'e_value': ".001",
                       'bitscore': "50",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "detect",
                       'vertical': "1",
                       'show_blanks': "1"
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
    
    
    ### Test 07: Specific Model with nucleotides (should fail)
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_07_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome_nuc_MSA()")
    def test_07_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome_nuc_MSA(self):
        test_name = 'test_07_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome_nuc_MSA'
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
                       'input_msa_refs': [self.MSA_refs[0]],  # Specific nuc MSA
                       'input_many_ref': self.genome_refs[0],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'sci_name',
                       'coalesce_output': 0,
                       'e_value': ".001",
                       'bitscore': "50",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "detect",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_Local_MSA_Group_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])
    
        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertEqual(report_obj['text_message'][0:7],"FAILURE")
        pass
    
    
    ### Test 08: Single Model against Single Genome, threshold above all hits
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_08_kb_hmmer_HMMER_MSA_Search_Genome_removeALL()")
    def test_08_kb_hmmer_HMMER_MSA_Search_Genome_removeALL(self):
        test_name = 'test_08_kb_hmmer_HMMER_MSA_Search_Genome_removeALL'
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
                       'input_msa_ref': self.MSA_refs[1],      # Single MSA
                       'input_many_ref': self.genome_refs[0],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_ver_sci_name',
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
        self.assertTrue(len(report_obj['objects_created']) == 0)
        pass
    
    
    ### Test 09: All Models in workspace against Single Genome, threshold above all hits
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_09_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome_removeALL()")
    def test_09_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome_removeALL(self):
        test_name = 'test_09_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome_removeALL'
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
                       'input_msa_refs': [self.MSA_refs[1], self.MSA_refs[2], self.MSA_refs[3]],  # Specific MSAs
                       'input_many_ref': self.genome_refs[0],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_ver_sci_name',
                       'coalesce_output': 0,
                       'e_value': ".001",
                       'bitscore': "500000000",
                       'overlap_fraction': "100.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "detect",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_Local_MSA_Group_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])
    
        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertTrue(len(report_obj['objects_created']) == 0)
        pass
    
    
    ### Test 09_02: All Models in workspace against Single Genome, threshold above all hits, show blanks
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_09_02_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome_removeALL_show_blanks()")
    def test_09_02_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome_removeALL_show_blanks(self):
        test_name = 'test_09_02_kb_hmmer_HMMER_Local_MSA_Group_Search_Genome_removeALL_show_blanks'
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
                       'input_msa_refs': [self.MSA_refs[1], self.MSA_refs[2], self.MSA_refs[3]],  # Specific MSAs
                       'input_many_ref': self.genome_refs[0],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_ver_sci_name',
                       'coalesce_output': 0,
                       'e_value': ".001",
                       'bitscore': "500000000",
                       'overlap_fraction': "100.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "detect",
                       'vertical': "1",
                       'show_blanks': "1"
                     }
        ret = self.getImpl().HMMER_Local_MSA_Group_Search(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['report_ref'])
    
        # check created obj
        #report_obj = self.getWsClient().get_objects2({'objects':[{'ref':ret['report_ref']}]})[0]['data']
        report_obj = self.getWsClient().get_objects([{'ref':ret['report_ref']}])[0]['data']
        self.assertTrue(len(report_obj['objects_created']) == 0)
        pass
    
    
    ### Test 10: dbCAN Models against Single Genome
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_10_kb_hmmer_HMMER_dbCAN_Search_Genome()")
    def test_10_kb_hmmer_HMMER_dbCAN_Search_Genome(self):
        test_name = 'test_10_kb_hmmer_HMMER_dbCAN_Search_Genome'
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
                       'input_dbCAN_AA_ids': ['ALL'],
                       'input_dbCAN_CBM_ids': ['ALL'],
                       'input_dbCAN_CE_ids': ['ALL'],
                       'input_dbCAN_GH_ids': ['ALL'],
                       'input_dbCAN_GT_ids': ['ALL'],
                       'input_dbCAN_PL_ids': ['ALL'],
                       'input_dbCAN_cellulosome_ids': ['ALL'],
                       'input_many_refs': [self.genome_refs[3]],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_ver_sci_name',
                       'coalesce_output': 0,
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "detect",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_dbCAN_Search(self.getContext(), parameters)[0]
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
    
    
    ### Test 10_02: dbCAN Models against Single AnnotatedMetagenomeAssembly
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_10_02_kb_hmmer_HMMER_dbCAN_Search_AMA()")
    def test_10_02_kb_hmmer_HMMER_dbCAN_Search_AMA(self):
        test_name = 'test_10_02_kb_hmmer_HMMER_dbCAN_Search_AMA'
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
                       'input_dbCAN_AA_ids': ['ALL'],
                       'input_dbCAN_CBM_ids': ['ALL'],
                       'input_dbCAN_CE_ids': ['ALL'],
                       'input_dbCAN_GH_ids': ['ALL'],
                       'input_dbCAN_GT_ids': ['ALL'],
                       'input_dbCAN_PL_ids': ['ALL'],
                       'input_dbCAN_cellulosome_ids': ['ALL'],
                       'input_many_refs': [self.ama_refs[0]],  # Single AnnotatedMetagenomeAssembly
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_ver_sci_name',
                       'coalesce_output': 0,
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       #'e_value': ".001",
                       'e_value': ".1",
                       #'bitscore': "50",
                       'bitscore': "5",
                       #'overlap_fraction': "50.0",
                       'overlap_fraction': "5.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_dbCAN_Search(self.getContext(), parameters)[0]
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
    
    
    ### Test 11: dbCAN Models against GenomeSet, DON'T coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_11_kb_hmmer_HMMER_dbCAN_Search_GenomeSet_NOcoalesce()")
    def test_11_kb_hmmer_HMMER_dbCAN_Search_GenomeSet_NOcoalesce(self):
        test_name = 'test_11_kb_hmmer_HMMER_dbCAN_Search_GenomeSet_NOcoalesce'
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
                       'input_dbCAN_AA_ids': ['AA1','AA6'],
                       'input_dbCAN_CBM_ids': ['CBM44','CBM48'],
                       'input_dbCAN_CE_ids': ['CE1'],
                       'input_dbCAN_GH_ids': ['GH5','GH13_20'],
                       'input_dbCAN_GT_ids': ['GT3','GT4'],
                       'input_dbCAN_PL_ids': ['PL9','PL22'],
                       'input_dbCAN_cellulosome_ids': ['dockerin','cohesin'],
                       'input_many_refs': [self.genomeSet_refs[0]],  # GenomeSet
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,  # KEY
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 0,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "0.1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_dbCAN_Search(self.getContext(), parameters)[0]
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
    
    
    ### Test 12: dbCAN Models against GenomeSet, DO coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_12_kb_hmmer_HMMER_dbCAN_Search_GenomeSet_coalesce()")
    def test_12_kb_hmmer_HMMER_dbCAN_Search_GenomeSet_coalesce(self):
        test_name = 'test_12_kb_hmmer_HMMER_dbCAN_Search_GenomeSet_coalesce'
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
                       'input_dbCAN_AA_ids': ['ALL'],
                       'input_dbCAN_CBM_ids': ['ALL'],
                       'input_dbCAN_CE_ids': ['ALL'],
                       'input_dbCAN_GH_ids': ['ALL'],
                       'input_dbCAN_GT_ids': ['ALL'],
                       'input_dbCAN_PL_ids': ['ALL'],
                       'input_dbCAN_cellulosome_ids': ['ALL'],
                       'input_many_refs': [self.genomeSet_refs[0]],  # GenomeSet
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 1,  # KEY
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'overlap_fraction': "50.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "detect",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_dbCAN_Search(self.getContext(), parameters)[0]
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

    
    ### Test 13: envbioelement Models against Single Genome
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_13_kb_hmmer_HMMER_envbioelement_Search_Genome()")
    def test_13_kb_hmmer_HMMER_envbioelement_Search_Genome(self):
        test_name = 'test_13_kb_hmmer_HMMER_envbioelement_Search_Genome'
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
                       'input_EnvBioelement_N_ids': ['ALL'],
                       'input_EnvBioelement_H_ids': ['ALL'],
                       'input_EnvBioelement_O_ids': ['ALL'],
                       'input_EnvBioelement_CFix_ids': ['ALL'],
                       'input_EnvBioelement_C1_ids': ['ALL'],
                       'input_EnvBioelement_CH4_ids': ['ALL'],
                       'input_EnvBioelement_CO_ids': ['ALL'],
                       'input_EnvBioelement_S_ids': ['ALL'],
                       'input_EnvBioelement_CN_ids': ['ALL'],
                       'input_EnvBioelement_CH4N2O_ids': ['ALL'],
                       'input_EnvBioelement_Se_ids': ['ALL'],
                       'input_EnvBioelement_Metal_ids': ['ALL'],
                       'input_EnvBioelement_As_ids': ['ALL'],
                       'input_EnvBioelement_Halo_ids': ['ALL'],
                       'input_many_refs': [self.genome_refs[3]],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'model_cov_perc': "35.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "detect",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_EnvBioelement_Search(self.getContext(), parameters)[0]
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


    ### Test 13_02: envbioelement Models against Single AMA
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_13_02_kb_hmmer_HMMER_envbioelement_Search_AMA()")
    def test_13_02_kb_hmmer_HMMER_envbioelement_Search_AMA(self):
        test_name = 'test_13_02_kb_hmmer_HMMER_envbioelement_Search_AMA'
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
                       'input_EnvBioelement_N_ids': ['ALL'],
                       'input_EnvBioelement_H_ids': ['ALL'],
                       'input_EnvBioelement_O_ids': ['ALL'],
                       'input_EnvBioelement_CFix_ids': ['ALL'],
                       'input_EnvBioelement_C1_ids': ['ALL'],
                       'input_EnvBioelement_CH4_ids': ['ALL'],
                       'input_EnvBioelement_CO_ids': ['ALL'],
                       'input_EnvBioelement_S_ids': ['ALL'],
                       'input_EnvBioelement_CN_ids': ['ALL'],
                       'input_EnvBioelement_CH4N2O_ids': ['ALL'],
                       'input_EnvBioelement_Se_ids': ['ALL'],
                       'input_EnvBioelement_Metal_ids': ['ALL'],
                       'input_EnvBioelement_As_ids': ['ALL'],
                       'input_EnvBioelement_Halo_ids': ['ALL'],
                       'input_many_refs': [self.ama_refs[0]],  # Single AMA
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'model_cov_perc': "35.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "detect",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_EnvBioelement_Search(self.getContext(), parameters)[0]
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


    ### Test 14: envbioelement Models against GenomeSet, DON'T coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_14_kb_hmmer_HMMER_envbioelement_Search_GenomeSet_NOcoalesce()")
    def test_14_kb_hmmer_HMMER_envbioelement_Search_GenomeSet_NOcoalesce(self):
        test_name = 'test_14_kb_hmmer_HMMER_envbioelement_Search_GenomeSet_NOcoalesce'
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
                       'input_EnvBioelement_N_ids': ['ALL'],
                       'input_EnvBioelement_H_ids': ['ALL'],
                       'input_EnvBioelement_O_ids': ['ALL'],
                       'input_EnvBioelement_CFix_ids': ['ALL'],
                       'input_EnvBioelement_C1_ids': ['ALL'],
                       'input_EnvBioelement_CH4_ids': ['ALL'],
                       'input_EnvBioelement_CO_ids': ['ALL'],
                       'input_EnvBioelement_S_ids': ['ALL'],
                       'input_EnvBioelement_CN_ids': ['ALL'],
                       'input_EnvBioelement_CH4N2O_ids': ['ALL'],
                       'input_EnvBioelement_Se_ids': ['ALL'],
                       'input_EnvBioelement_Metal_ids': ['ALL'],
                       'input_EnvBioelement_As_ids': ['ALL'],
                       'input_EnvBioelement_Halo_ids': ['ALL'],
                       'input_many_refs': [self.genomeSet_refs[0]],  # GenomeSet
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,  # KEY
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'model_cov_perc': "35.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_EnvBioelement_Search(self.getContext(), parameters)[0]
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


    ### Test 14_02: envbioelement Models against GenomeSet+AMA, DON'T coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_14_02_kb_hmmer_HMMER_envbioelement_Search_GenomeSet_AMA_NOcoalesce()")
    def test_14_02_kb_hmmer_HMMER_envbioelement_Search_GenomeSet_AMA_NOcoalesce(self):
        test_name = 'test_14_02_kb_hmmer_HMMER_envbioelement_Search_GenomeSet_AMA_NOcoalesce'
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
                       'input_EnvBioelement_N_ids': ['nirB'],
                       'input_EnvBioelement_H_ids': ['NONE'],
                       'input_EnvBioelement_O_ids': ['NONE'],
                       'input_EnvBioelement_CFix_ids': ['NONE'],
                       'input_EnvBioelement_C1_ids': ['NONE'],
                       'input_EnvBioelement_CH4_ids': ['NONE'],
                       'input_EnvBioelement_CO_ids': ['NONE'],
                       'input_EnvBioelement_S_ids': ['NONE'],
                       'input_EnvBioelement_CN_ids': ['NONE'],
                       'input_EnvBioelement_CH4N2O_ids': ['NONE'],
                       'input_EnvBioelement_Se_ids': ['NONE'],
                       'input_EnvBioelement_Metal_ids': ['NONE'],
                       'input_EnvBioelement_As_ids': ['NONE'],
                       'input_EnvBioelement_Halo_ids': ['NONE'],
                       'input_many_refs': [self.genomeSet_refs[0], self.ama_refs[0]], # GenomeSet+AMA
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,  # KEY
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'model_cov_perc': "35.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_EnvBioelement_Search(self.getContext(), parameters)[0]
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


    ### Test 15: envbioelement Models against GenomeSet, DO coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_15_kb_hmmer_HMMER_envbioelement_Search_GenomeSet_coalesce()")
    def test_15_kb_hmmer_HMMER_envbioelement_Search_GenomeSet_coalesce(self):
        test_name = 'test_15_kb_hmmer_HMMER_envbioelement_Search_GenomeSet_coalesce'
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
                       'input_EnvBioelement_N_ids': ['ALL'],
                       'input_EnvBioelement_H_ids': ['ALL'],
                       'input_EnvBioelement_O_ids': ['ALL'],
                       'input_EnvBioelement_CFix_ids': ['ALL'],
                       'input_EnvBioelement_C1_ids': ['ALL'],
                       'input_EnvBioelement_CH4_ids': ['ALL'],
                       'input_EnvBioelement_CO_ids': ['ALL'],
                       'input_EnvBioelement_S_ids': ['ALL'],
                       'input_EnvBioelement_CN_ids': ['ALL'],
                       'input_EnvBioelement_CH4N2O_ids': ['ALL'],
                       'input_EnvBioelement_Se_ids': ['ALL'],
                       'input_EnvBioelement_Metal_ids': ['ALL'],
                       'input_EnvBioelement_As_ids': ['ALL'],
                       'input_EnvBioelement_Halo_ids': ['ALL'],
                       'input_many_refs': [self.genomeSet_refs[0]],  # GenomeSet
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_ver_sci_name',
                       'coalesce_output': 1,  # KEY
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'model_cov_perc': "35.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "0.1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_EnvBioelement_Search(self.getContext(), parameters)[0]
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


    ### Test 16: PhyloMarkers Models against Single Genome
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_16_kb_hmmer_HMMER_PhyloMarkers_Search_Genome()")
    def test_16_kb_hmmer_HMMER_PhyloMarkers_Search_Genome(self):
        test_name = 'test_16_kb_hmmer_HMMER_PhyloMarkers_Search_Genome'
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
                       'input_PhyloMarkers_Univ_ids': ['ALL'],
                       'input_PhyloMarkers_B_ribo_pol_ids': ['ALL'],
                       'input_PhyloMarkers_B_other_ids': ['ALL'],
                       'input_PhyloMarkers_A_ribo_pol_ids': ['ALL'],
                       'input_PhyloMarkers_A_other_ids': ['ALL'],
                       'input_many_refs': [self.genome_refs[3]],  # Single Genome
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'model_cov_perc': "35.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "detect",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_PhyloMarkers_Search(self.getContext(), parameters)[0]
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


    ### Test 16_02: PhyloMarkers Models against Single AMA
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_16_02_kb_hmmer_HMMER_PhyloMarkers_Search_AMA()")
    def test_16_02_kb_hmmer_HMMER_PhyloMarkers_Search_AMA(self):
        test_name = 'test_16_02_kb_hmmer_HMMER_PhyloMarkers_Search_AMA'
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
                       'input_PhyloMarkers_Univ_ids': ['ALL'],
                       'input_PhyloMarkers_B_ribo_pol_ids': ['ALL'],
                       'input_PhyloMarkers_B_other_ids': ['ALL'],
                       'input_PhyloMarkers_A_ribo_pol_ids': ['ALL'],
                       'input_PhyloMarkers_A_other_ids': ['ALL'],
                       'input_many_refs': [self.ama_refs[0]],  # Single AMA
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'model_cov_perc': "35.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "detect",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_PhyloMarkers_Search(self.getContext(), parameters)[0]
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


    ### Test 17: PhyloMarkers Models against GenomeSet, DON'T coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_17_kb_hmmer_HMMER_PhyloMarkers_Search_GenomeSet_NOcoalesce()")
    def test_17_kb_hmmer_HMMER_PhyloMarkers_Search_GenomeSet_NOcoalesce(self):
        test_name = 'test_17_kb_hmmer_HMMER_PhyloMarkers_Search_GenomeSet_NOcoalesce'
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
                       'input_PhyloMarkers_Univ_ids': ['ALL'],
                       'input_PhyloMarkers_B_ribo_pol_ids': ['ALL'],
                       'input_PhyloMarkers_B_other_ids': ['ALL'],
                       'input_PhyloMarkers_A_ribo_pol_ids': ['ALL'],
                       'input_PhyloMarkers_A_other_ids': ['ALL'],
                       'input_many_refs': [self.genomeSet_refs[0]],  # GenomeSet
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,  # KEY
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'model_cov_perc': "35.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_PhyloMarkers_Search(self.getContext(), parameters)[0]
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


    ### Test 17_02: PhyloMarkers Models against GenomeSet+AMA, DON'T coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_17_02_kb_hmmer_HMMER_PhyloMarkers_Search_GenomeSet_AMA_NOcoalesce()")
    def test_17_02_kb_hmmer_HMMER_PhyloMarkers_Search_GenomeSet_AMA_NOcoalesce(self):
        test_name = 'test_17_02_kb_hmmer_HMMER_PhyloMarkers_Search_GenomeSet_AMA_NOcoalesce'
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
                       'input_PhyloMarkers_Univ_ids': ['ALL'],
                       'input_PhyloMarkers_B_ribo_pol_ids': ['ALL'],
                       'input_PhyloMarkers_B_other_ids': ['ALL'],
                       'input_PhyloMarkers_A_ribo_pol_ids': ['ALL'],
                       'input_PhyloMarkers_A_other_ids': ['ALL'],
                       'input_many_refs': [self.genomeSet_refs[0], self.ama_refs[0]], # GenomeSet+AMA
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,  # KEY
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'model_cov_perc': "35.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_PhyloMarkers_Search(self.getContext(), parameters)[0]
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


    ### Test 17_03: PhyloMarkers Models against AMA+AMA, DON'T coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_17_03_kb_hmmer_HMMER_PhyloMarkers_Search_AMA_AMA_NOcoalesce()")
    def test_17_03_kb_hmmer_HMMER_PhyloMarkers_Search_AMA_AMA_NOcoalesce(self):
        test_name = 'test_17_03_kb_hmmer_HMMER_PhyloMarkers_Search_AMA_AMA_NOcoalesce'
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
                       'input_PhyloMarkers_Univ_ids': ['ALL'],
                       'input_PhyloMarkers_B_ribo_pol_ids': ['ALL'],
                       'input_PhyloMarkers_B_other_ids': ['ALL'],
                       'input_PhyloMarkers_A_ribo_pol_ids': ['ALL'],
                       'input_PhyloMarkers_A_other_ids': ['ALL'],
                       'input_many_refs': [self.ama_refs[0], self.ama_refs[1]], # AMA+AMA
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_sci_name',
                       'coalesce_output': 0,  # KEY
                       'show_target_block_headers': 0,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'model_cov_perc': "35.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_PhyloMarkers_Search(self.getContext(), parameters)[0]
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


    ### Test 18: PhyloMarkers Models against GenomeSet, DO coalesce output
    #
    # uncomment to skip this test
    # HIDE @unittest.skip("skipped test test_18_kb_hmmer_HMMER_PhyloMarkers_Search_GenomeSet_coalesce()")
    def test_18_kb_hmmer_HMMER_PhyloMarkers_Search_GenomeSet_coalesce(self):
        test_name = 'test_18_kb_hmmer_HMMER_PhyloMarkers_Search_GenomeSet_coalesce'
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
                       'input_PhyloMarkers_Univ_ids': ['ALL'],
                       'input_PhyloMarkers_B_ribo_pol_ids': ['ALL'],
                       'input_PhyloMarkers_B_other_ids': ['ALL'],
                       'input_PhyloMarkers_A_ribo_pol_ids': ['ALL'],
                       'input_PhyloMarkers_A_other_ids': ['ALL'],
                       'input_many_refs': [self.genomeSet_refs[0]],  # GenomeSet
                       'output_filtered_name': obj_out_name,
                       'genome_disp_name_config': 'obj_name_ver_sci_name',
                       'coalesce_output': 1,  # KEY
                       'show_target_block_headers': 1,
                       'save_ALL_featureSets': 1,
                       'save_ANY_featureSets': 1,
                       'e_value': ".001",
                       'bitscore': "50",
                       'model_cov_perc': "35.0",
                       'maxaccepts': "1000",
                       'heatmap': "1",
                       'low_val': "0.1",
                       'vertical': "1",
                       'show_blanks': "0"
                     }
        ret = self.getImpl().HMMER_PhyloMarkers_Search(self.getContext(), parameters)[0]
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
