package uk.ac.ebi.pride.tools.pride_spectra_clustering.consensus_spectrum_builder;

import java.util.List;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

/**
 * Creates a consensus spectrum from
 * the passed list of spectra.
 * @author jg
 *
 */
public interface ConsensusSpectrumBuilder {
	/**
	 * Creates a consensus spectrum based on the
	 * passed list of spectra.
	 * @param spectra A list of spectra as sorted peak lists according to intensity.
	 * @return A list of Peaks sorted according to their intensities.
	 */
	public List<Peak> buildConsensusSpectrum(List<List<Peak>> spectra);
}
