package us.kbase.kbhmmer;

import com.fasterxml.jackson.core.type.TypeReference;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import us.kbase.auth.AuthToken;
import us.kbase.common.service.JsonClientCaller;
import us.kbase.common.service.JsonClientException;
import us.kbase.common.service.RpcContext;
import us.kbase.common.service.UnauthorizedException;

/**
 * <p>Original spec-file module name: kb_hmmer</p>
 * <pre>
 * ** A KBase module: kb_hmmer
 * **
 * ** This module contains HMMER Hidden Markov Model Sequence Search and Alignment
 * **
 * </pre>
 */
public class KbHmmerClient {
    private JsonClientCaller caller;
    private String serviceVersion = null;


    /** Constructs a client with a custom URL and no user credentials.
     * @param url the URL of the service.
     */
    public KbHmmerClient(URL url) {
        caller = new JsonClientCaller(url);
    }
    /** Constructs a client with a custom URL.
     * @param url the URL of the service.
     * @param token the user's authorization token.
     * @throws UnauthorizedException if the token is not valid.
     * @throws IOException if an IOException occurs when checking the token's
     * validity.
     */
    public KbHmmerClient(URL url, AuthToken token) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, token);
    }

    /** Constructs a client with a custom URL.
     * @param url the URL of the service.
     * @param user the user name.
     * @param password the password for the user name.
     * @throws UnauthorizedException if the credentials are not valid.
     * @throws IOException if an IOException occurs when checking the user's
     * credentials.
     */
    public KbHmmerClient(URL url, String user, String password) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, user, password);
    }

    /** Constructs a client with a custom URL
     * and a custom authorization service URL.
     * @param url the URL of the service.
     * @param user the user name.
     * @param password the password for the user name.
     * @param auth the URL of the authorization server.
     * @throws UnauthorizedException if the credentials are not valid.
     * @throws IOException if an IOException occurs when checking the user's
     * credentials.
     */
    public KbHmmerClient(URL url, String user, String password, URL auth) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, user, password, auth);
    }

    /** Get the token this client uses to communicate with the server.
     * @return the authorization token.
     */
    public AuthToken getToken() {
        return caller.getToken();
    }

    /** Get the URL of the service with which this client communicates.
     * @return the service URL.
     */
    public URL getURL() {
        return caller.getURL();
    }

    /** Set the timeout between establishing a connection to a server and
     * receiving a response. A value of zero or null implies no timeout.
     * @param milliseconds the milliseconds to wait before timing out when
     * attempting to read from a server.
     */
    public void setConnectionReadTimeOut(Integer milliseconds) {
        this.caller.setConnectionReadTimeOut(milliseconds);
    }

    /** Check if this client allows insecure http (vs https) connections.
     * @return true if insecure connections are allowed.
     */
    public boolean isInsecureHttpConnectionAllowed() {
        return caller.isInsecureHttpConnectionAllowed();
    }

    /** Deprecated. Use isInsecureHttpConnectionAllowed().
     * @deprecated
     */
    public boolean isAuthAllowedForHttp() {
        return caller.isAuthAllowedForHttp();
    }

    /** Set whether insecure http (vs https) connections should be allowed by
     * this client.
     * @param allowed true to allow insecure connections. Default false
     */
    public void setIsInsecureHttpConnectionAllowed(boolean allowed) {
        caller.setInsecureHttpConnectionAllowed(allowed);
    }

    /** Deprecated. Use setIsInsecureHttpConnectionAllowed().
     * @deprecated
     */
    public void setAuthAllowedForHttp(boolean isAuthAllowedForHttp) {
        caller.setAuthAllowedForHttp(isAuthAllowedForHttp);
    }

    /** Set whether all SSL certificates, including self-signed certificates,
     * should be trusted.
     * @param trustAll true to trust all certificates. Default false.
     */
    public void setAllSSLCertificatesTrusted(final boolean trustAll) {
        caller.setAllSSLCertificatesTrusted(trustAll);
    }
    
    /** Check if this client trusts all SSL certificates, including
     * self-signed certificates.
     * @return true if all certificates are trusted.
     */
    public boolean isAllSSLCertificatesTrusted() {
        return caller.isAllSSLCertificatesTrusted();
    }
    /** Sets streaming mode on. In this case, the data will be streamed to
     * the server in chunks as it is read from disk rather than buffered in
     * memory. Many servers are not compatible with this feature.
     * @param streamRequest true to set streaming mode on, false otherwise.
     */
    public void setStreamingModeOn(boolean streamRequest) {
        caller.setStreamingModeOn(streamRequest);
    }

    /** Returns true if streaming mode is on.
     * @return true if streaming mode is on.
     */
    public boolean isStreamingModeOn() {
        return caller.isStreamingModeOn();
    }

    public void _setFileForNextRpcResponse(File f) {
        caller.setFileForNextRpcResponse(f);
    }

    public String getServiceVersion() {
        return this.serviceVersion;
    }

    public void setServiceVersion(String newValue) {
        this.serviceVersion = newValue;
    }

    /**
     * <p>Original spec-file function name: HMMER_MSA_Search</p>
     * <pre>
     * Method for HMMER search of an MSA against many sequences 
     * **
     * **    overloading as follows:
     * **        input_msa_ref: MSA
     * **        input_many_ref: SequenceSet, FeatureSet, Genome, GenomeSet, AMA (note: SequenceSet deactivated)
     * **        output_name: SequenceSet (if input_many is SequenceSet), (else) FeatureSet
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbhmmer.HMMERParams HMMERParams} (original type "HMMER_Params")
     * @return   instance of type {@link us.kbase.kbhmmer.HMMEROutput HMMEROutput} (original type "HMMER_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public HMMEROutput hMMERMSASearch(HMMERParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<HMMEROutput>> retType = new TypeReference<List<HMMEROutput>>() {};
        List<HMMEROutput> res = caller.jsonrpcCall("kb_hmmer.HMMER_MSA_Search", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: HMMER_Local_MSA_Group_Search</p>
     * <pre>
     * Method for HMMER search of a Local MSA Group (found automatically within workspace) against many sequences 
     * **
     * **    overloading as follows:
     * **        input_many_ref: SequenceSet, FeatureSet, Genome, GenomeSet, AMA (note: SeqeuenceSet deactivated)
     * **        output_name: SequenceSet (if input_many is SequenceSet), (else) FeatureSet
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbhmmer.HMMERLocalMSAGroupParams HMMERLocalMSAGroupParams} (original type "HMMER_Local_MSA_Group_Params")
     * @return   instance of type {@link us.kbase.kbhmmer.HMMEROutput HMMEROutput} (original type "HMMER_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public HMMEROutput hMMERLocalMSAGroupSearch(HMMERLocalMSAGroupParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<HMMEROutput>> retType = new TypeReference<List<HMMEROutput>>() {};
        List<HMMEROutput> res = caller.jsonrpcCall("kb_hmmer.HMMER_Local_MSA_Group_Search", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: HMMER_dbCAN_Search</p>
     * <pre>
     * Method for HMMER search of dbCAN Markov Models of CAZy families
     * **
     * **    overloading as follows:
     * **        input_many_ref: SequenceSet, FeatureSet, Genome, GenomeSet, AMA (note: SequenceSet deactivated)
     * **        output_name: SequenceSet (if input_many is SequenceSet), (else) FeatureSet
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbhmmer.HMMERDbCANParams HMMERDbCANParams} (original type "HMMER_dbCAN_Params")
     * @return   instance of type {@link us.kbase.kbhmmer.HMMEROutput HMMEROutput} (original type "HMMER_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public HMMEROutput hMMERDbCANSearch(HMMERDbCANParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<HMMEROutput>> retType = new TypeReference<List<HMMEROutput>>() {};
        List<HMMEROutput> res = caller.jsonrpcCall("kb_hmmer.HMMER_dbCAN_Search", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: HMMER_EnvBioelement_Search</p>
     * <pre>
     * Method for HMMER search of Markov Models of environmental bioelement families
     * **
     * **    overloading as follows:
     * **        input_many_ref: SequenceSet, FeatureSet, Genome, GenomeSet, AMA (note: SeqquenceSet deactivated)
     * **        output_name: SequenceSet (if input_many is SequenceSet), (else) FeatureSet
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbhmmer.HMMEREnvBioelementParams HMMEREnvBioelementParams} (original type "HMMER_EnvBioelement_Params")
     * @return   instance of type {@link us.kbase.kbhmmer.HMMEROutput HMMEROutput} (original type "HMMER_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public HMMEROutput hMMEREnvBioelementSearch(HMMEREnvBioelementParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<HMMEROutput>> retType = new TypeReference<List<HMMEROutput>>() {};
        List<HMMEROutput> res = caller.jsonrpcCall("kb_hmmer.HMMER_EnvBioelement_Search", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: HMMER_PhyloMarkers_Search</p>
     * <pre>
     * Method for HMMER search of Markov Models of phylogenetic marker families
     * **
     * **    overloading as follows:
     * **        input_many_ref: SequenceSet, FeatureSet, Genome, GenomeSet, AMA (note: SeqquenceSet deactivated)
     * **        output_name: SequenceSet (if input_many is SequenceSet), (else) FeatureSet
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbhmmer.HMMERPhyloMarkersParams HMMERPhyloMarkersParams} (original type "HMMER_PhyloMarkers_Params")
     * @return   instance of type {@link us.kbase.kbhmmer.HMMEROutput HMMEROutput} (original type "HMMER_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public HMMEROutput hMMERPhyloMarkersSearch(HMMERPhyloMarkersParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<HMMEROutput>> retType = new TypeReference<List<HMMEROutput>>() {};
        List<HMMEROutput> res = caller.jsonrpcCall("kb_hmmer.HMMER_PhyloMarkers_Search", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    public Map<String, Object> status(RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        TypeReference<List<Map<String, Object>>> retType = new TypeReference<List<Map<String, Object>>>() {};
        List<Map<String, Object>> res = caller.jsonrpcCall("kb_hmmer.status", args, retType, true, false, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }
}
