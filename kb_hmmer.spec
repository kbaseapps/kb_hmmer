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

	float e_value;
	float bitscore;
	float maxaccepts;

/*	float ident_thresh;
	float overlap_fraction;
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
	

    /*  Methods for HMMER search of an MSA against many sequences 
    **
    **    overloading as follows:
    **        input_msa_ref: MSA
    **        input_many_ref: SingleEndLibrary, FeatureSet, Genome, GenomeSet
    **        output_name: SingleEndLibrary (if input_many is SELib), (else) FeatureSet
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

	bool  coalesce_output;
	float e_value;
	float bitscore;
	float maxaccepts;
/*	float ident_thresh;
	float overlap_fraction;
*/
        bool  heatmap;
	bool  vertical;  /* only supports true for now */
	bool  show_blanks;
    } HMMER_Local_MSA_Group_Params;


    /*  Methods for HMMER search of a Local MSA Group (found automatically within workspace) against many sequences 
    **
    **    overloading as follows:
    **        input_many_ref: SingleEndLibrary, FeatureSet, Genome, GenomeSet
    **        output_name: SingleEndLibrary (if input_many is SELib), (else) FeatureSet
    */
    funcdef HMMER_Local_MSA_Group_Search (HMMER_Local_MSA_Group_Params params)  returns (HMMER_Output) authentication required;
};
