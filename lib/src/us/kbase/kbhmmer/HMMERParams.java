
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
 * <p>Original spec-file type: HMMER_Params</p>
 * <pre>
 * HMMER Input Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "input_many_ref",
    "input_msa_ref",
    "output_filtered_name",
    "e_value",
    "bitscore",
    "overlap_perc",
    "maxaccepts"
})
public class HMMERParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("input_many_ref")
    private String inputManyRef;
    @JsonProperty("input_msa_ref")
    private String inputMsaRef;
    @JsonProperty("output_filtered_name")
    private String outputFilteredName;
    @JsonProperty("e_value")
    private Double eValue;
    @JsonProperty("bitscore")
    private Double bitscore;
    @JsonProperty("overlap_perc")
    private Double overlapPerc;
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

    public HMMERParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
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

    public HMMERParams withInputManyRef(String inputManyRef) {
        this.inputManyRef = inputManyRef;
        return this;
    }

    @JsonProperty("input_msa_ref")
    public String getInputMsaRef() {
        return inputMsaRef;
    }

    @JsonProperty("input_msa_ref")
    public void setInputMsaRef(String inputMsaRef) {
        this.inputMsaRef = inputMsaRef;
    }

    public HMMERParams withInputMsaRef(String inputMsaRef) {
        this.inputMsaRef = inputMsaRef;
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

    public HMMERParams withOutputFilteredName(String outputFilteredName) {
        this.outputFilteredName = outputFilteredName;
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

    public HMMERParams withEValue(Double eValue) {
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

    public HMMERParams withBitscore(Double bitscore) {
        this.bitscore = bitscore;
        return this;
    }

    @JsonProperty("overlap_perc")
    public Double getOverlapPerc() {
        return overlapPerc;
    }

    @JsonProperty("overlap_perc")
    public void setOverlapPerc(Double overlapPerc) {
        this.overlapPerc = overlapPerc;
    }

    public HMMERParams withOverlapPerc(Double overlapPerc) {
        this.overlapPerc = overlapPerc;
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

    public HMMERParams withMaxaccepts(Double maxaccepts) {
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
        return ((((((((((((((((((("HMMERParams"+" [workspaceName=")+ workspaceName)+", inputManyRef=")+ inputManyRef)+", inputMsaRef=")+ inputMsaRef)+", outputFilteredName=")+ outputFilteredName)+", eValue=")+ eValue)+", bitscore=")+ bitscore)+", overlapPerc=")+ overlapPerc)+", maxaccepts=")+ maxaccepts)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
