
package us.kbase.kbhmmer;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: HMMER_Local_MSA_Group_Params</p>
 * <pre>
 * HMMER Local MSA Group Input Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "use_all_local_MSAs",
    "input_msa_refs",
    "input_many_ref",
    "output_filtered_name",
    "coalesce_output",
    "e_value",
    "bitscore",
    "maxaccepts"
})
public class HMMERLocalMSAGroupParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("use_all_local_MSAs")
    private Long useAllLocalMSAs;
    @JsonProperty("input_msa_refs")
    private String inputMsaRefs;
    @JsonProperty("input_many_ref")
    private String inputManyRef;
    @JsonProperty("output_filtered_name")
    private String outputFilteredName;
    @JsonProperty("coalesce_output")
    private Long coalesceOutput;
    @JsonProperty("e_value")
    private Double eValue;
    @JsonProperty("bitscore")
    private Double bitscore;
    @JsonProperty("maxaccepts")
    private Double maxaccepts;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public HMMERLocalMSAGroupParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("use_all_local_MSAs")
    public Long getUseAllLocalMSAs() {
        return useAllLocalMSAs;
    }

    @JsonProperty("use_all_local_MSAs")
    public void setUseAllLocalMSAs(Long useAllLocalMSAs) {
        this.useAllLocalMSAs = useAllLocalMSAs;
    }

    public HMMERLocalMSAGroupParams withUseAllLocalMSAs(Long useAllLocalMSAs) {
        this.useAllLocalMSAs = useAllLocalMSAs;
        return this;
    }

    @JsonProperty("input_msa_refs")
    public String getInputMsaRefs() {
        return inputMsaRefs;
    }

    @JsonProperty("input_msa_refs")
    public void setInputMsaRefs(String inputMsaRefs) {
        this.inputMsaRefs = inputMsaRefs;
    }

    public HMMERLocalMSAGroupParams withInputMsaRefs(String inputMsaRefs) {
        this.inputMsaRefs = inputMsaRefs;
        return this;
    }

    @JsonProperty("input_many_ref")
    public String getInputManyRef() {
        return inputManyRef;
    }

    @JsonProperty("input_many_ref")
    public void setInputManyRef(String inputManyRef) {
        this.inputManyRef = inputManyRef;
    }

    public HMMERLocalMSAGroupParams withInputManyRef(String inputManyRef) {
        this.inputManyRef = inputManyRef;
        return this;
    }

    @JsonProperty("output_filtered_name")
    public String getOutputFilteredName() {
        return outputFilteredName;
    }

    @JsonProperty("output_filtered_name")
    public void setOutputFilteredName(String outputFilteredName) {
        this.outputFilteredName = outputFilteredName;
    }

    public HMMERLocalMSAGroupParams withOutputFilteredName(String outputFilteredName) {
        this.outputFilteredName = outputFilteredName;
        return this;
    }

    @JsonProperty("coalesce_output")
    public Long getCoalesceOutput() {
        return coalesceOutput;
    }

    @JsonProperty("coalesce_output")
    public void setCoalesceOutput(Long coalesceOutput) {
        this.coalesceOutput = coalesceOutput;
    }

    public HMMERLocalMSAGroupParams withCoalesceOutput(Long coalesceOutput) {
        this.coalesceOutput = coalesceOutput;
        return this;
    }

    @JsonProperty("e_value")
    public Double getEValue() {
        return eValue;
    }

    @JsonProperty("e_value")
    public void setEValue(Double eValue) {
        this.eValue = eValue;
    }

    public HMMERLocalMSAGroupParams withEValue(Double eValue) {
        this.eValue = eValue;
        return this;
    }

    @JsonProperty("bitscore")
    public Double getBitscore() {
        return bitscore;
    }

    @JsonProperty("bitscore")
    public void setBitscore(Double bitscore) {
        this.bitscore = bitscore;
    }

    public HMMERLocalMSAGroupParams withBitscore(Double bitscore) {
        this.bitscore = bitscore;
        return this;
    }

    @JsonProperty("maxaccepts")
    public Double getMaxaccepts() {
        return maxaccepts;
    }

    @JsonProperty("maxaccepts")
    public void setMaxaccepts(Double maxaccepts) {
        this.maxaccepts = maxaccepts;
    }

    public HMMERLocalMSAGroupParams withMaxaccepts(Double maxaccepts) {
        this.maxaccepts = maxaccepts;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((((((((("HMMERLocalMSAGroupParams"+" [workspaceName=")+ workspaceName)+", useAllLocalMSAs=")+ useAllLocalMSAs)+", inputMsaRefs=")+ inputMsaRefs)+", inputManyRef=")+ inputManyRef)+", outputFilteredName=")+ outputFilteredName)+", coalesceOutput=")+ coalesceOutput)+", eValue=")+ eValue)+", bitscore=")+ bitscore)+", maxaccepts=")+ maxaccepts)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
