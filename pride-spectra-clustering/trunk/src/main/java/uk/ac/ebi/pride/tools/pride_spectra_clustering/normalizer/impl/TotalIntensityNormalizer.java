package uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.impl;

import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.IntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

/**
 * Normalizes a spectrum's intensities so that
 * the spectrum's total intensity matches a given
 * number.
 * 
 * This method was used by Frank et al. (2008) JPR: Clustering millions of spectra
 * To calculate the dot-product they then used
 * 1+ln(intensity) as the peak's intensity
 * @author jg
 *
 */
public class TotalIntensityNormalizer implements IntensityNormalizer {
	private double totalIntensity = 1000;

	public List<Peak> normalizeSpectrum(List<Peak> spectrum) {
		// get the max intensity
		Double specTotalIntensity = 0.0;
		
		for (Peak p : spectrum)
			specTotalIntensity += p.getIntensity();
		
		// if there's no suitable max intensity, return the unchanged spectrum
		if (specTotalIntensity <= 0)
			return spectrum;
		
		// calculate the ratio
		double ratio = totalIntensity / specTotalIntensity;
		
		// create the new spectrum
		List<Peak> normalizedSpectrum = new ArrayList<Peak>(spectrum.size());
		
		for (Peak p : spectrum)
			normalizedSpectrum.add(new Peak(p.getMz(), p.getIntensity() * ratio));
		
		return normalizedSpectrum;
	}

	/**
	 * Set the total intensity the spectrum's
	 * intensity should be set to.
	 * @param totalIntensity
	 */
	public void setTotalIntensity(double totalIntensity) {
		this.totalIntensity = totalIntensity;
	}
}
