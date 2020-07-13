/*
** A KBase module: kb_hmmer
**
** This module contains HMMER Hidden Markov Model Sequence Search and Alignment
**
*/

module kb_hmmer {

    /*
    ** The workspace object refs are of form:
    **
    **    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_name']}])
    **
    ** "ref" means the entire name combining the workspace id and the object name
    ** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
    ** "name" is a string identifier of a workspace or object.  This is received from Narrative.
    */
    typedef string workspace_name;
    typedef string sequence;
    typedef string data_obj_name;
    typedef string data_obj_ref;
    typedef int    bool;


    /* HMMER Input Params
    */
    typedef structure {
        workspace_name workspace_name;
/*	sequence       input_one_sequence;
	data_obj_ref   input_one_ref;
*/
	data_obj_ref   input_many_ref;
	data_obj_ref   input_msa_ref;
        data_obj_name  output_filtered_name;
        string         genome_disp_name_config;

	float e_value;
	float bitscore;
	float overlap_perc;
	float maxaccepts;

/*	float ident_thresh;
*/
    } HMMER_Params;


    /* HMMER Output
    */
    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
/*       data_obj_ref  output_filtered_ref;
*
*        int n_initial_seqs;
*        int n_seqs_matched;
*        int n_seqs_notmatched;
*/
    } HMMER_Output;
	

    /*  Method for HMMER search of an MSA against many sequences 
    **
    **    overloading as follows:
    **        input_msa_ref: MSA
    **        input_many_ref: SequenceSet, FeatureSet, Genome, GenomeSet
    **        output_name: SequenceSet (if input_many is SequenceSet), (else) FeatureSet
    */
    funcdef HMMER_MSA_Search (HMMER_Params params)  returns (HMMER_Output) authentication required;


    /* HMMER Local MSA Group Input Params
    */
    typedef structure {
        workspace_name workspace_name;
/*	sequence       input_one_sequence;
	data_obj_ref   input_one_ref;
*/
	data_obj_ref   input_msa_refs;
	data_obj_ref   input_many_ref;
        data_obj_name  output_filtered_name;
        string         genome_disp_name_config;

	bool  coalesce_output;
	float e_value;
	float bitscore;
	float overlap_perc;
	float maxaccepts;
/*	float ident_thresh;
*/
        bool  heatmap;
	bool  vertical;  /* only supports true for now */
	bool  show_blanks;
    } HMMER_Local_MSA_Group_Params;


    /*  Method for HMMER search of a Local MSA Group (found automatically within workspace) against many sequences 
    **
    **    overloading as follows:
    **        input_many_ref: SequenceSet, FeatureSet, Genome, GenomeSet
    **        output_name: SequenceSet (if input_many is SequenceSet), (else) FeatureSet
    */
    funcdef HMMER_Local_MSA_Group_Search (HMMER_Local_MSA_Group_Params params)  returns (HMMER_Output) authentication required;


    /* HMMER dbCAN Input Params
    */
    typedef structure {
        workspace_name workspace_name;
/*	sequence       input_one_sequence;
	data_obj_ref   input_one_ref;
*/
	data_obj_ref   input_dbCAN_AA_ids;
	data_obj_ref   input_dbCAN_CBM_ids;
	data_obj_ref   input_dbCAN_CE_ids;
	data_obj_ref   input_dbCAN_GH_ids;
	data_obj_ref   input_dbCAN_GT_ids;
	data_obj_ref   input_dbCAN_PL_ids;
	data_obj_ref   input_dbCAN_cellulosome_ids;
	data_obj_ref   input_many_ref;
        data_obj_name  output_filtered_name;
        string         genome_disp_name_config;

	bool  coalesce_output;
	float e_value;
	float bitscore;
	float overlap_perc;
	float maxaccepts;
/*	float ident_thresh;
*/
        bool  heatmap;
	bool  vertical;  /* only supports true for now */
	bool  show_blanks;
    } HMMER_dbCAN_Params;


    /*  Method for HMMER search of dbCAN Markov Models of CAZy families
    **
    **    overloading as follows:
    **        input_many_ref: SequenceSet, FeatureSet, Genome, GenomeSet
    **        output_name: SequenceSet (if input_many is SequenceSet), (else) FeatureSet
    */
    funcdef HMMER_dbCAN_Search (HMMER_dbCAN_Params params)  returns (HMMER_Output) authentication required;


    /* HMMER EnvBioelement Input Params
    */
    typedef structure {
        workspace_name workspace_name;
/*  sequence       input_one_sequence;
    data_obj_ref   input_one_ref;
*/
    data_obj_ref   input_env-bioelement_N_ids;
    data_obj_ref   input_env-bioelement_S_ids;
    data_obj_ref   input_env-bioelement_O_ids;
    data_obj_ref   input_env-bioelement_CH4_ids;
    data_obj_ref   input_env-bioelement_CFix_ids;
    data_obj_ref   input_env-bioelement_CMono_ids;
    data_obj_ref   input_env-bioelement_C1_ids;
    data_obj_ref   input_env-bioelement_H_ids;
    data_obj_ref   input_env-bioelement_Halo_ids;
    data_obj_ref   input_env-bioelement_As_ids;
    data_obj_ref   input_env-bioelement_Se_ids;
    data_obj_ref   input_env-bioelement_Ur_ids;
    data_obj_ref   input_env-bioelement_Me_ids;
    data_obj_ref   input_env-bioelement_CN_ids;
    data_obj_ref   input_many_ref;
        data_obj_name  output_filtered_name;
        string         genome_disp_name_config;

    bool  coalesce_output;
    float e_value;
    float bitscore;
    float overlap_perc;
    float maxaccepts;
/*  float ident_thresh;
*/
        bool  heatmap;
    bool  vertical;  /* only supports true for now */
    bool  show_blanks;
    } HMMER_EnvBioelement_Params;


    /*  Method for HMMER search of Markov Models of environmental bioelement families
    **
    **    overloading as follows:
    **        input_many_ref: SequenceSet, FeatureSet, Genome, GenomeSet
    **        output_name: SequenceSet (if input_many is SequenceSet), (else) FeatureSet
    */
    funcdef HMMER_EnvBioelement_Search (HMMER_EnvBioelement_Params params)  returns (HMMER_Output) authentication required;
};
