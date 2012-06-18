package uk.ac.ebi.pride.tools.cluster.pride_spectra_clustering_example;

import java.io.File;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.SpectraClustering;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.impl.FrankEtAlClustering;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.ClusteringSpectrum;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.SpectraCluster;

/**
 * Hello world!
 *
 */
public class App 
{
	private final static Pattern peptideSequencePattern = Pattern.compile(".*sequence=([A-Z-]+).*");
	
    public static void main( String[] args )
    {
        System.out.println("+-------------------------------------+");
        System.out.println("|   PRIDE Spectra Clustering example  |");
        System.out.println("+-------------------------------------+\n");
        
        try {
        	// use the first argument as filename
        	if (args.length < 1) {
        		System.out.println("Usage: java -jar pride-spectra-clustering-example-1.0-SNAPSHOT.jar [MGF file]");
        		System.exit(1);
        		return;
        	}
        	
        	String filename = args[0];
        	
        	DecimalFormatSymbols decimalSymbols = new DecimalFormatSymbols(Locale.US);
        	DecimalFormat decimalFormat = new DecimalFormat("#.##", decimalSymbols);
        	
        	// load the spectra from the MGF file
        	System.out.print("Loading MGF file...");
        	
        	// WARNING: this does not work if everything is packed
        	//MgfFile mgfFile = new MgfFile(new File(spectraFile.toURI()));
        	MgfFile mgfFile = new MgfFile(new File(filename));
        	// load all spectra into memory
        	List<ClusteringSpectrum> spectra = loadSpectraFromFile(mgfFile);
        	
        	System.out.println("Done. (" + mgfFile.getMs2QueryCount() + " spectra loaded).");
        	
        	// Initialize the clustering
        	SpectraClustering clustering = new FrankEtAlClustering();
        	// set the variables required for the clustering process
        	// the settings used here are the same ones as the ones
        	// used to cluster the PRIDE database
        	clustering.setClusteringRounds(4);
        	clustering.setSimilarityThreshold(0.7);
        	
        	// cluster the spectra
        	System.out.print("Clustering spectra...");
        	long start = System.currentTimeMillis();
        	List<SpectraCluster> createdCluster = clustering.clusterConvertedSpectra(spectra);
        	long duration = System.currentTimeMillis() - start;
        	
        	// show the results
        	System.out.println("Done (" + decimalFormat.format(duration / 1000) + " seconds, " + createdCluster.size() + " cluster created.)\n\n");
        	
        	System.out.println("--------------- RESULTS -----------------");

        	for (SpectraCluster cluster : createdCluster) {
        		System.out.println("\n- Cluster (" + decimalFormat.format(cluster.getAverageMz()) + " m/z) -");
        		
        		// count the sequences
        		Map<String, Integer> sequences = new HashMap<String, Integer>();
        		
        		// iterate over all spectra and extract the peptide's sequence
        		// from the spectrum's id
        		for (ClusteringSpectrum spectrum : cluster.getSpectra()) {
        			Matcher sequenceMatcher = peptideSequencePattern.matcher(spectrum.getId());
        			
        			if (!sequenceMatcher.find())
        				throw new Exception("Failed to extract peptide sequence from " + spectrum.getId());
        			
        			String sequence = sequenceMatcher.group(1);
        			
        			if (!sequences.containsKey(sequence))
        				sequences.put(sequence, 0);
        			
        			sequences.put(sequence, sequences.get(sequence) + 1);
        		}
        		
        		// print the ratios of the sequences in the cluster
        		for (String sequence : sequences.keySet()) {
        			System.out.println(sequence + " (" + sequences.get(sequence) + " / " + cluster.getSpectra().size() + " spectra)");
        		}
        	}
        }
        catch (Exception e) {
        	System.out.println("Failed.");
        	e.printStackTrace();
        	System.exit(1);
        }
    }
    
    /**
	 * Loads all spectra from the MGF file and returns
	 * them as a List of Spectra. The spectra's id is
	 * replaced by their title.
	 * @param mgfFile
	 * @return
	 * @throws Exception 
	 */
	private static List<ClusteringSpectrum> loadSpectraFromFile(MgfFile mgfFile) throws Exception {
		List<ClusteringSpectrum> spectra = new ArrayList<ClusteringSpectrum>(mgfFile.getMs2QueryCount());
		
		Iterator<Ms2Query> it = mgfFile.getMs2QueryIterator();
		Set<String> processedIds = new HashSet<String>();
		
		while(it.hasNext()) {
			Ms2Query query = it.next();
			
			// make sure every spectrum is only used once
			if (processedIds.contains(query.getTitle()))
				continue;
			
			processedIds.add(query.getTitle());
			
			// set the intensity to 1 in case it's missing
			if (query.getPeptideIntensity() == null)
				query.setPeptideIntensity(1.0);
			
			// calculate the charge in case it's missing
			if (query.getPrecursorCharge() == null) {
				throw new Exception("Spectrum is missing precursor charge.");
			}
			
			// change the id to title
			ClusteringSpectrum spectrum = new ClusteringSpectrum(
					query.getTitle(), 
					query.getPrecursorMZ(), 
					query.getPrecursorIntensity(), 
					query.getPrecursorCharge(), 
					query.getPeakList(), 
					null, 
					query.getMsLevel());
			
			spectra.add(spectrum);
		}
		
		return spectra;
	}
}
