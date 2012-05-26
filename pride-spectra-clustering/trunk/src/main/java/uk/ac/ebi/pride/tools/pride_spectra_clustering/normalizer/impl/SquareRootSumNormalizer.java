package uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.impl;

import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.IntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

/**
 * Normalizes a spectrum by dividing
 * the intensities through the square
 * root of the sum of all squared
 * intensities.
 * @author jg
 *
 */
public class SquareRootSumNormalizer implements IntensityNormalizer {

	public List<Peak> normalizeSpectrum(List<Peak> spectrum) {
		// get the intensities
		List<Double> intensities = new ArrayList<Double>(spectrum.size());
		for (Peak p : spectrum)
			intensities.add(p.getIntensity());
		
		// calculate the sum of the squared intensities
		double sumSquaredIntensities = 0;
		
		for (Double intensity : intensities)
			sumSquaredIntensities += Math.pow(intensity, 2);
		
		// calculate the ratio by taking the square root of the summed intensities
		double ratio = Math.sqrt(sumSquaredIntensities);
		
		// create the new spectrum
		List<Peak> normalizedSpectrum = new ArrayList<Peak>(spectrum.size());
		
		for (Peak p : spectrum)
			normalizedSpectrum.add(new Peak(p.getMz(), p.getIntensity() * ratio));
		
		return normalizedSpectrum;
	}

}
