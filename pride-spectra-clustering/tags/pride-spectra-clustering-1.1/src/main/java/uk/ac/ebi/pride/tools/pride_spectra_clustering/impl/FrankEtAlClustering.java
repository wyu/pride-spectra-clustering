package uk.ac.ebi.pride.tools.pride_spectra_clustering.impl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.log4j.Logger;

import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.SpectraClustering;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.consensus_spectrum_builder.ConsensusSpectrumBuilder;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.consensus_spectrum_builder.impl.FrankEtAlConsensusSpectrumBuilder;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.IntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.impl.TotalIntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.quality_check.QualityChecker;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.quality_check.impl.SignalToNoiseChecker;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.similarity_checker.SimilarityChecker;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.similarity_checker.impl.FrankEtAlDotProduct;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.ClusteringSpectrum;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.SpectraCluster;

public class FrankEtAlClustering implements SpectraClustering {
	/**
	 * Indicates whether info about the processing times
	 * should be printed to STDOUT.
	 */
	public static boolean PRINT_STATS_INFO = false;
	/**
	 * Logger to use
	 */
	private final static Logger logger = Logger.getLogger(FrankEtAlClustering.class);
	/**
	 * Maximum number of clustering rounds
	 */
	private int maxClusteringRounds = 10;
	/**
	 * Similarity checker to use.
	 */
	public static final SimilarityChecker similarityChecker = new FrankEtAlDotProduct();
	/**
	 * The similarity threshold to use to cluster spectra
	 */
	private double similarityThreshold = 0.6;
	/**
	 * The range to use when precursor m/z values
	 * are considered equal.
	 */
	private final double precursorRange = 2;
	/**
	 * The quality checker is only used to sort the spectra.
	 */
	private final QualityChecker qualityChecker = new SignalToNoiseChecker();
	/**
	 * The consensus spectrum builder to use.
	 */
	public static final ConsensusSpectrumBuilder consensusBuilder = new FrankEtAlConsensusSpectrumBuilder();
	/**
	 * The intensity normalizer used
	 */
	public static final IntensityNormalizer normalizer = new TotalIntensityNormalizer();

	/**
	 * Cluster the passed List of Spectra
	 * and returns the resulting clusters 
	 * as a List of SpectraClusterS.
	 * @param spectra List of spectra to cluster
	 * @return A List of resulting SpectraCluster
	 */
	public List<SpectraCluster> clusterSpectra(List<Spectrum> inputSpectra) {
		// convert all spectra into sorted peak lists
		List<ClusteringSpectrum> spectra = new ArrayList<ClusteringSpectrum>(inputSpectra.size());
		for (Spectrum s : inputSpectra)
			spectra.add(new ClusteringSpectrum(s));
		
		return clusterConvertedSpectra(spectra);
	}
	
	@Override
	public List<SpectraCluster> clusterConvertedSpectra(
			List<ClusteringSpectrum> spectra) {		
		return clusterConvertedSpectra(spectra, new ArrayList<SpectraCluster>());
	}

	@Override
	public List<SpectraCluster> clusterConvertedSpectra(
			List<ClusteringSpectrum> spectra,
			List<SpectraCluster> cluster) {
		// sort the spectra
		Collections.sort(spectra, new SpectraComparator());
		
		// cluster the spectra
		long start = System.currentTimeMillis();
		clusterSpectra(cluster, spectra);
		long duration = System.currentTimeMillis() - start;
		if (PRINT_STATS_INFO) System.out.println("Initial clustering: " + duration + " msec");
		
		int currentRound = 0;
		
		List<ClusteringSpectrum> nonFittingSpectra = null;
		
		while ((nonFittingSpectra == null || nonFittingSpectra.size() > 0) && currentRound <= maxClusteringRounds) {
			// merge similar cluster
			start = System.currentTimeMillis();
			int clusterBefore = cluster.size();
			cluster = mergeCluster(cluster);
			duration = System.currentTimeMillis() - start;
			if (PRINT_STATS_INFO) System.out.println("Merging " + currentRound + ": " + duration + " msec (" + (clusterBefore - cluster.size()) + " merged)");
			
			// remove non-fitting peptides from the cluster
			start = System.currentTimeMillis();
			nonFittingSpectra = removeNonFittingSpectra(cluster);
			duration = System.currentTimeMillis() - start;
			if (PRINT_STATS_INFO) System.out.println("Non fitting " + currentRound + ": " + duration + " msec (" + nonFittingSpectra.size() + ")");
			
			// cluster the non-fitting peptides again
			start = System.currentTimeMillis();
			clusterSpectra(cluster, nonFittingSpectra);
			duration = System.currentTimeMillis() - start;
			if (PRINT_STATS_INFO) System.out.println("Clustering " + currentRound + ": " + duration + " msec");
			
			currentRound++;
		}
		
		// return only cluster containing spectra
		List<SpectraCluster> nonEmptyCluster = new ArrayList<SpectraCluster>(cluster.size());
		
		for (SpectraCluster c : cluster) {
			if (c != null && c.getClusterSize() > 0)
				nonEmptyCluster.add(c);
		}
		
		return nonEmptyCluster;
	}

	/**
	 * Removes all spectra from the cluster
	 * that no longer fit the consensus spectrum
	 * and returns them as a list.
	 * @param cluster
	 * @return
	 */
	private List<ClusteringSpectrum> removeNonFittingSpectra(List<SpectraCluster> cluster) {
		List<ClusteringSpectrum> nonFittingSpectra = new ArrayList<ClusteringSpectrum>();
		
		for (SpectraCluster c : cluster) {
			List<ClusteringSpectrum> removedSpectra = c.removeNonFittingSpectra(similarityChecker, similarityThreshold);
			
			if (removedSpectra.size() > 0) {
				logger.debug(c.getAverageMz() + " contained " + removedSpectra.size() + " non-fitting spectra");
				nonFittingSpectra.addAll(removedSpectra);
			}
		}
		
		return nonFittingSpectra;
	}

	/**
	 * Clusters the passed spectra with the
	 * existing clusters or in case they don't
	 * fit any of the existing ones, creates a
	 * new cluster.
	 * @param cluster
	 * @param spectra
	 */
	private void clusterSpectra(List<SpectraCluster> cluster, List<ClusteringSpectrum> spectra) {
		// loop through all spectra
		for (ClusteringSpectrum s : spectra) {			
			SpectraCluster mostSimilarCluster = null;
			double highestSimilarity = 0;
			
			// check if there's a cluster where the spectrum fits
			for (SpectraCluster c : cluster) {
				double similarity = similarityChecker.assessSimilarity(c.getConsensusSpectrum(), s.getPeaklist(), c.getAverageMz(), s.getPrecursorMZ(), c.getAverageCharge(), s.getPrecursorCharge().doubleValue());
				
				if (similarity >= similarityThreshold && similarity > highestSimilarity) {
					highestSimilarity = similarity;
					mostSimilarCluster = c;
				}
			}
			
			// check if the spectrum fits a cluster
			if (mostSimilarCluster != null) {
				mostSimilarCluster.addSpectrum(s);
			}
			// if the spectrum wasn't added to any cluster, create a new one
			else {
				SpectraCluster newCluster = new SpectraCluster(consensusBuilder, normalizer);
				newCluster.addSpectrum(s);
				cluster.add(newCluster);
			}
		}
	}

	/**
	 * Merges all cluster all cluster based
	 * on whether their consensus spectra
	 * are similar and returns a list of
	 * merged cluster.
	 * WARNING: The passed list is altered!
	 * @param cluster
	 * @return
	 */
	private List<SpectraCluster> mergeCluster(List<SpectraCluster> cluster) {
		// now try to cluster the cluster with each other - until no more cluster can be merged
		boolean clusterMerged = true;
		
		while (clusterMerged) {
			clusterMerged = false; 
			
			for (int i = 0; i < cluster.size(); i++) {
				SpectraCluster c1 = cluster.get(i);
				
				// ignore deleted cluster
				if (c1 == null)
					continue;
				
				for (int j = 0; j < cluster.size(); j++) {
					// ignore the same cluster
					if (j == i)
						continue;
					
					SpectraCluster c2 = cluster.get(j);
					
					// ignore deleted cluster
					if (c2 == null)
						continue;
					
					double similarity = similarityChecker.assessSimilarity(c1.getConsensusSpectrum(), c2.getConsensusSpectrum(), c1.getAverageMz(), c2.getAverageMz(), c1.getAverageCharge(), c2.getAverageCharge());
					
					if (similarity >= similarityThreshold) {
						logger.debug("Merging cluster " + c1.getAverageMz() + " with " + c2.getAverageMz()); 
						clusterMerged = true;
						
						// add the spectra from one cluster to the other
						c1.addSpectra(c2.getSpectra());
						
						// delete the second cluster
						cluster.set(j, null);
					}
				}
			}
		}
			
		// return a list without any null or empty cluster
		List<SpectraCluster> cleanList = new ArrayList<SpectraCluster>(cluster.size());
		
		for (SpectraCluster c : cluster) {
			if (c != null && c.getClusterSize() > 0)
				cleanList.add(c);
		}
		
		return cleanList;
	}

	/**
	 * Used to sort spectra first according
	 * to the precursor m/z and then according
	 * to their "quality".
	 * @author jg
	 *
	 */
	private class SpectraComparator implements Comparator<ClusteringSpectrum> {

		public int compare(ClusteringSpectrum spec1, ClusteringSpectrum spec2) {
			// first check whether the precursor m/z is different
			if (spec1.getPrecursorMZ() - precursorRange > spec2.getPrecursorMZ())
				return -1;
			if (spec1.getPrecursorMZ() + precursorRange < spec2.getPrecursorMZ())
				return 1;
			
			// as the m/z ranges are the same check the quality
			double quality1 = qualityChecker.assessQuality(spec1.getPeaklist());
			double quality2 = qualityChecker.assessQuality(spec2.getPeaklist());
			
			if (quality1 > quality2)
				return -1;
			if (quality1 < quality2)
				return 1;
			
			return 0;
		}
		
	}

	public Double getSimilarityThreshold() {
		return similarityThreshold;
	}

	public void setSimilarityThreshold(Double similarityThreshold) {
		this.similarityThreshold = similarityThreshold;
	}

	public int getClusteringRounds() {
		return maxClusteringRounds;
	}

	public void setClusteringRounds(int rounds) {
		maxClusteringRounds = rounds;
	}

	public String getDescription() {
		String description = 
			"FrankEtAlClustering: adapted version " +
			"of the clustering algorithm described by Frank et al. in Nat. Meth. (2011).\n" +
			"Threshold: " + similarityThreshold + "\n" +
			"Clustering Rounds: " + maxClusteringRounds;
		
		return description;
	}
}
