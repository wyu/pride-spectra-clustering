package uk.ac.ebi.pride.tools.pride_spectra_clustering.quality_check.impl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.quality_check.QualityChecker;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

/**
 * Assesses a spectrum's signal-to-noise
 * ratio as a measure of the spectrum's
 * quality. The method used is taken from
 * Lam et al. (2008), Nat Methods 5(10):873-875 
 * @author jg
 *
 */
public class SignalToNoiseChecker implements QualityChecker {

	/**
	 * Calculates a spectrum's signal-to-noise ratio
	 * by taking the 2nd-6th highest peak's intensity
	 * and dividing it by the median intensity of all
	 * peaks.
	 * This method should only be used on normalized
	 * spectra.
	 */
	public double assessQuality(List<Peak> spectrum) {
		// get the intensities
		ArrayList<Double> intensities = new ArrayList<Double>(spectrum.size());
		for (Peak p : spectrum)
			intensities.add(p.getIntensity());
		
		Collections.sort(intensities);
		
		// get the total intensity of the 2nd-6th highest peak
		double highestPeakIntensity = 0.0;
		
		// TODO: react to small spectra
		if (intensities.size() < 6)
			return 0.0;
		
		for (int n = intensities.size() - 6; n < intensities.size() - 1; n++)
			highestPeakIntensity += intensities.get(n);
		
		highestPeakIntensity = highestPeakIntensity / 5;
		
		// get the median
		double nPeaks = intensities.size();
		double median = 0.0;
		
		// check if there's an even number of peaks
		if (nPeaks %2 == 0)
			median = (intensities.get((int) nPeaks / 2 - 1) + intensities.get((int) nPeaks / 2)) / 2;
		else
			median = intensities.get((int) (nPeaks / 2));
		
		return highestPeakIntensity / median;
	}

	public boolean requiresNormalization() {
		return true;
	}
}
