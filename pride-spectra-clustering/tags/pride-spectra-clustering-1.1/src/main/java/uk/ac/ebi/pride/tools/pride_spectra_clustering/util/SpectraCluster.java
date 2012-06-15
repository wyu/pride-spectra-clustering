package uk.ac.ebi.pride.tools.pride_spectra_clustering.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.log4j.Logger;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.consensus_spectrum_builder.ConsensusSpectrumBuilder;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.IntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.similarity_checker.SimilarityChecker;

/**
 * Represents a cluster of spectra
 * 
 * @author jg
 * 
 */
public class SpectraCluster {
	/**
	 * The logger to use.
	 */
	private final static Logger logger = Logger.getLogger(SpectraCluster.class);
	/**
	 * The consensus spectrum builder used to create the consensus spectrum.
	 */
	private final ConsensusSpectrumBuilder consensusSpectrumBuilder;
	/**
	 * The intensity normalizer to use to normalize the single spectrum's
	 * intensities.
	 */
	private final IntensityNormalizer normalizer;
	/**
	 * The cluster's consensus spectrum.
	 */
	protected List<Peak> consensusSpectrum;
	/**
	 * A list of all the spectra that were added to the cluster.
	 */
	protected List<ClusteringSpectrum> addedSpectra;
	/**
	 * The average m/z of this cluster.
	 */
	protected double averageMz = 0;
	/**
	 * The cluster's average charge.
	 */
	protected double averageCharge = 0;
	/**
	 * The total precursor intensity.
	 */
	protected double totalIntensity = 0;
	/**
	 * Additional information about the cluster
	 */
	protected String additional;

	/**
	 * Creates a new SpectraCluster.
	 * 
	 * @param consensusSpectrumBuilder
	 *            The ConensusSpectrumBuilder used for this cluster.
	 * @param normalizer
	 *            The IntensityNormalizer to use for the clustering.
	 */
	public SpectraCluster(ConsensusSpectrumBuilder consensusSpectrumBuilder,
			IntensityNormalizer normalizer) {
		this.consensusSpectrumBuilder = consensusSpectrumBuilder;
		this.normalizer = normalizer;
	}

	/**
	 * Initialize a SpectraCluster object that was loaded from a previous
	 * clustering process.
	 * 
	 * @param consensusSpectrumBuilder
	 * @param normalizer
	 * @param consensusSpectrum
	 * @param addedSpectra
	 * @param averageMz
	 * @param averageCharge
	 * @param totalIntensity
	 */
	public SpectraCluster(ConsensusSpectrumBuilder consensusSpectrumBuilder, 
			IntensityNormalizer normalizer, List<Peak> consensusSpectrum,
			List<ClusteringSpectrum> addedSpectra, 
			double averageMz,
			double averageCharge, double totalIntensity) {
		this(consensusSpectrumBuilder, normalizer);

		// rebuild the consensus spectrum
		this.consensusSpectrum = consensusSpectrum;		
		this.addedSpectra = addedSpectra;
		this.averageMz = averageMz;
		this.averageCharge = averageCharge;
		this.totalIntensity = totalIntensity;
	}

	/**
	 * Returns the consensus spectrum representing the cluster of spectra.
	 * 
	 * @return A Map representing the consensus spectrum of this cluster with
	 *         the m/z values as keys and their intensities as values
	 */
	public List<Peak> getConsensusSpectrum() {
		return consensusSpectrum;
	}
	
	/**
	 * Add multiple spectra to the cluster at once.
	 * @param spectra
	 */
	public void addSpectra(List<ClusteringSpectrum> spectra) throws IllegalArgumentException {
		// make sure the spectra are valid
		for (ClusteringSpectrum spectrum : spectra) {
			// make sure the spectrum contains all the required fields.
			if (spectrum.getPrecursorCharge() == null)
				throw new IllegalArgumentException(
						"The spectrum's charge must not be null.");
			if (spectrum.getPrecursorIntensity() == null)
				throw new IllegalArgumentException(
						"The spectrum's intensity must not be null.");
			if (spectrum.getPrecursorMZ() == null)
				throw new IllegalArgumentException(
						"The spectrum's precursor m/z must not be null.");
		}
		
		// make sure the addedSpectra list was initialized
		if (addedSpectra == null)
			addedSpectra = new ArrayList<ClusteringSpectrum>();
		
		double weightedMz = 0, totalIntensity = 0, totalCharge = 0;
		
		for (ClusteringSpectrum spectrum : spectra) {
			// calculate the weighted mz
			if (weightedMz == 0)
				weightedMz = spectrum.getPrecursorMZ();
			else {
				weightedMz = 
					( (weightedMz * totalIntensity) + (spectrum.getPrecursorMZ() * spectrum.getPrecursorIntensity()) )
					/ (totalIntensity + spectrum.getPrecursorIntensity());
			}
			
			// update the average intensity
			totalIntensity += spectrum.getPrecursorIntensity();
			// update the average charge
			totalCharge += spectrum.getPrecursorCharge().doubleValue();
		}
		
		// update the cluster's charge, weighted m/z, and intensity
		this.averageCharge = ((this.averageCharge * addedSpectra.size()) + (totalCharge)) / (addedSpectra.size() + spectra.size());
		this.averageMz = ((weightedMz * totalIntensity) + (this.averageMz * this.totalIntensity)) / (totalIntensity + this.totalIntensity);
		this.totalIntensity += totalIntensity;
		
		// add the spectra
		addedSpectra.addAll(spectra);
		
		buildConsensusSpectrum();
	}

	/**
	 * Adds the passed spectrum to the cluster. This function does not check any
	 * similarities or similar.
	 * 
	 * @param spectrum
	 */
	public void addSpectrum(ClusteringSpectrum spectrum)
			throws IllegalArgumentException {
		// make sure the spectrum contains all the required fields.
		if (spectrum.getPrecursorCharge() == null)
			throw new IllegalArgumentException(
					"The spectrum's charge must not be null.");
		if (spectrum.getPrecursorIntensity() == null)
			throw new IllegalArgumentException(
					"The spectrum's intensity must not be null.");
		if (spectrum.getPrecursorMZ() == null)
			throw new IllegalArgumentException(
					"The spectrum's precursor m/z must not be null.");

		// make sure the addedSpectra list was initialized
		if (addedSpectra == null)
			addedSpectra = new ArrayList<ClusteringSpectrum>();

		// add the spectrum
		addedSpectra.add(spectrum);

		// update the average mz using the weighted precursor intensity
		double weightedMzCluster = averageMz * totalIntensity
				/ (totalIntensity + spectrum.getPrecursorIntensity());
		double weightedMzSpectrum = spectrum.getPrecursorMZ()
				* spectrum.getPrecursorIntensity()
				/ (totalIntensity + spectrum.getPrecursorIntensity());
		averageMz = weightedMzCluster + weightedMzSpectrum;
		// update the average intensity
		totalIntensity += spectrum.getPrecursorIntensity();
		// update the average charge
		averageCharge = averageCharge * (addedSpectra.size() - 1)
				/ addedSpectra.size()
				+ spectrum.getPrecursorCharge().doubleValue()
				/ addedSpectra.size();

		buildConsensusSpectrum();
	}

	private void buildConsensusSpectrum() {		
		// rebuild the consensus spectrum
		List<List<Peak>> peakLists = new ArrayList<List<Peak>>(
				addedSpectra.size());

		// create the list of normalized peak lists
		for (ClusteringSpectrum s : addedSpectra)
			peakLists.add(s.getPeaklist());

		consensusSpectrum = normalizer
				.normalizeSpectrum(consensusSpectrumBuilder
						.buildConsensusSpectrum(peakLists));
	}

	/**
	 * Removes the given spectrum from the cluster.
	 * 
	 * @param spectrum
	 */
	public void removeSpectrum(ClusteringSpectrum spectrum) {
		if (addedSpectra == null)
			return;

		addedSpectra.remove(spectrum);

		// rebuild the consensus spectrum
		buildConsensusSpectrum();
	}

	public double getAverageMz() {
		return averageMz;
	}

	public double getAverageCharge() {
		return averageCharge;
	}

	public double getAverageIntensity() {
		return totalIntensity;
	}

	public int getClusterSize() {
		if (addedSpectra == null)
			return 0;

		return addedSpectra.size();
	}

	public List<ClusteringSpectrum> getSpectra() {
		return addedSpectra;
	}

	@Override
	public String toString() {
		if (addedSpectra == null)
			return "[]";

		String string = "[";
		for (ClusteringSpectrum s : addedSpectra)
			string += (string.length() > 1 ? "," : "") + s.getId();

		string += "]";

		return string;
	}

	/**
	 * Removes any spectra from the cluster that are no longer fitting the
	 * consensus spectrum.
	 * 
	 * @return The List of Spectra that are no longer fitting the cluster's
	 *         consensus spectrum
	 */
	public List<ClusteringSpectrum> removeNonFittingSpectra(
			SimilarityChecker checker, double similariyThreshold) {
		List<ClusteringSpectrum> removedSpectra = new ArrayList<ClusteringSpectrum>();

		for (int i = 0; i < addedSpectra.size(); i++) {
			ClusteringSpectrum s = addedSpectra.get(i);

			double similarity = checker.assessSimilarity(consensusSpectrum, s
					.getPeaklist(), averageMz, s.getPrecursorMZ(),
					averageCharge, (s.getPrecursorCharge() != null) ? s
							.getPrecursorCharge().doubleValue() : null);

			if (similarity < similariyThreshold) {
				// add the spectrum to the list of removed spectra
				removedSpectra.add(s);
				// set the spectrum to null in the array
				addedSpectra.set(i, null);
			}
		}

		// make sure spectra were removed
		if (removedSpectra.size() < 1)
			return removedSpectra;

		// "clean" the addedSpectra by removing all null values
		List<ClusteringSpectrum> cleanAddedSpectra = new ArrayList<ClusteringSpectrum>(
				addedSpectra.size() - removedSpectra.size());

		for (ClusteringSpectrum s : addedSpectra) {
			if (s != null)
				cleanAddedSpectra.add(s);
		}

		// replace the added spectra with the "cleaned" version
		addedSpectra = cleanAddedSpectra;

		// make sure there are spectra left
		if (addedSpectra.size() < 1) {
			averageCharge = 0;
			averageMz = 0;
			consensusSpectrum = Collections.emptyList();
			return removedSpectra;
		}

		// update the average charge and mz (weighted by intensity)
		double totalIntensities = 0, mzTimesIntensities = 0, totalCharge = 0;

		for (ClusteringSpectrum s : addedSpectra) {
			totalIntensities += s.getPrecursorIntensity();
			mzTimesIntensities += s.getPrecursorMZ()
					* s.getPrecursorIntensity();
			totalCharge += s.getPrecursorCharge();
		}

		averageMz = mzTimesIntensities / totalIntensities;
		averageCharge = totalCharge / addedSpectra.size();

		// rebuild the consensus spectrum
		buildConsensusSpectrum();

		logger.debug("Removed " + removedSpectra.size()
				+ " non-fitting spectra. averageMz = " + averageMz
				+ ", averageCharge = " + averageCharge);

		return removedSpectra;
	}

	public String getAdditional() {
		return additional;
	}

	public void setAdditional(String additional) {
		this.additional = additional;
	}
}
