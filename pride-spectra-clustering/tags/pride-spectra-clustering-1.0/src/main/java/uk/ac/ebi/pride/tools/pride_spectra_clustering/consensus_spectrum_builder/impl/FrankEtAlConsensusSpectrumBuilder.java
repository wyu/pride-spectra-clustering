package uk.ac.ebi.pride.tools.pride_spectra_clustering.consensus_spectrum_builder.impl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.consensus_spectrum_builder.ConsensusSpectrumBuilder;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.PeakIntensityComparator;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.PeakMzComparator;

/**
 * Generates a consensus spectrum as described
 * by Frank etl al. (2008) JPR.
 * @author jg
 *
 */
public class FrankEtAlConsensusSpectrumBuilder implements
		ConsensusSpectrumBuilder {
	/**
	 * The final m/z threshold to use to combine
	 * peaks as identical.
	 */
	private double finalMzThreshold = 0.4;
	/**
	 * The step size to use when iteratively
	 * merging identical peaks. This is not
	 * done with the finalMzThreshold from the
	 * beginning but starts from mzThresholdStep
	 * towards finalMzThreshold.
	 */
	private double mzThresholdStep = 0.1;

	public List<Peak> buildConsensusSpectrum(List<List<Peak>> spectra) {
		// make sure something sensible was passed
		if (spectra == null || spectra.size() < 1)
			return null;
		
		// if there's only one spectrum in the list, return this spectrum
		if (spectra.size() == 1)
			return spectra.get(0);
		
		// initialize the consensus spectrum - expect 2000 peaks for the beginning
		Map<Double, Peak> consensusSpectrum = new HashMap<Double, Peak>(2000);
		
		// add the peaks from all spectra to the consensus spectrum
		addAllPeaks(consensusSpectrum, spectra);
		
		// merge identical peaks in the consensus spectrum
		List<Peak> mergedConsensusSpectrum = mergeIdenticalPeaks(consensusSpectrum);
		
		// adapt the peak intensities using the following formula: I = I * (0.95 + 0.05 * (1 + pi))^5 where pi is the peaks probability
		mergedConsensusSpectrum = adaptPeakIntensities(mergedConsensusSpectrum, spectra.size());
		
		// filter the spectrum
		List<Peak> filteredSpectrum = filterSpectrum(mergedConsensusSpectrum);
		
		// sort the spectrum according to intensities
		Collections.sort(filteredSpectrum, PeakIntensityComparator.getInstance());
		
		return filteredSpectrum;
	}

	/**
	 * Filters the passed spectrum keeping only the
	 * top 5 peaks per 100 Da
	 * @param mergedConsensusSpectrum
	 * @return
	 */
	private List<Peak> filterSpectrum(
			List<Peak> mergedConsensusSpectrum) {
		// expect to keep 1% - just a wild guess
		List<Peak> filteredSpectrum = new ArrayList<Peak>(mergedConsensusSpectrum.size() / 100);
		
		// sort the passed spectrum
		Collections.sort(mergedConsensusSpectrum, PeakMzComparator.getInstance());
		
		// process the peaks using the sliding window		
		for (double startMz = 0, endMz = 100; endMz <= 5000; endMz += 100, startMz += 100) {
			List<Peak> peakBuffer = new ArrayList<Peak>();
			
			// fill the peak buffer with all peaks within that range
			for (Peak p : mergedConsensusSpectrum) {
				if (p.getMz() < startMz)
					continue;
				if (p.getMz() > endMz)
					break;
				
				peakBuffer.add(p);
			}
			
			// sort the buffer
			Collections.sort(peakBuffer, PeakIntensityComparator.getInstance());
			
			// take the 5 highest peaks
			for (int i = peakBuffer.size() - 1, counter = 0; i >= 0 && counter < 5; i--, counter++)
				filteredSpectrum.add(peakBuffer.get(i));
		}
		
		return filteredSpectrum;
	}

	private List<Peak> adaptPeakIntensities(List<Peak> mergedConsensusSpectrum, double numberOfSpectra) {
		List<Peak> adaptedSpectrum = new ArrayList<Peak>(mergedConsensusSpectrum.size());
		
		for (Peak p : mergedConsensusSpectrum) {
			if (p == null)
				continue;
			
			double peakProbability = (double) p.getCount() / numberOfSpectra;
			double newIntensity = p.getIntensity() * (0.95 + 0.05 * Math.pow(1 + peakProbability,5) );
			
			adaptedSpectrum.add(new Peak(p.getMz(), newIntensity, p.getCount()));
		}
		
		return adaptedSpectrum;
	}

	/**
	 * Adds all peaks from the passed List of
	 * spectra to the consensus spectrum.
	 * @param consensusSpectrum
	 * @param spectra
	 */
	private void addAllPeaks(Map<Double, Peak> consensusSpectrum,
			List<List<Peak>> spectra) {
		// make sure the parameters are valid
		if (spectra == null || spectra.size() < 1)
			return;
		if (consensusSpectrum == null)
			return;
		
		// process the spectra
		for (List<Peak> spectrum : spectra) {
			// process the current spectrum's peaks
			for (Peak peak : spectrum) {
				// ignore 0 intensity peaks
				if (peak.getIntensity() == 0)
					continue;
				
				// if the consensus spectrum doesn't contain the peak yet, simply add it
				if (!consensusSpectrum.containsKey(peak.getMz())) {
					consensusSpectrum.put(peak.getMz(), peak);
				}
				else {
					// if the peak already exists, sum it up
					consensusSpectrum.put(peak.getMz(), new Peak(
							peak.getMz(), 
							consensusSpectrum.get(peak.getMz()).getIntensity() + peak.getIntensity(), 
							consensusSpectrum.get(peak.getMz()).getCount() + 1)
					);
				}
			}
		}
	}
	
	/**
	 * 
	 * @param consensusSpectrum
	 */
	private List<Peak> mergeIdenticalPeaks(Map<Double, Peak> consensusSpectrum) {
		// convert the spectrum into a list of Peaks
		List<Peak> peaks = new ArrayList<Peak>( consensusSpectrum.values() );
		
		// based on the set range
		for (double range = mzThresholdStep; range <= finalMzThreshold; range += mzThresholdStep) {
			// sort the list according to m/z values
			Collections.sort(peaks, PeakMzComparator.getInstance());
			
			// as the list is sorted, peaks only have to be checked in one "direction"
			for (int i = 0; i < peaks.size() - 1; i++) {
				Peak current = peaks.get(i);
				Peak next    = peaks.get(i + 1);
				
				if (current == null || next == null)
					continue;
				
				// check if the next peak falls within the range
				if (next.getMz() <= current.getMz() + range) {
					// calculate the new weighted m/z
					double weightedMz = (next.getIntensity() * next.getMz() + current.getIntensity() * current.getMz()) / (next.getIntensity() + current.getIntensity());
					
					Peak mergedPeak = new Peak(weightedMz, current.getIntensity() + next.getIntensity(), current.getCount() + next.getCount());
					
					// remove the current peak from the array
					peaks.set(i, null);
					// set the next peak to the merged one
					peaks.set(i + 1, mergedPeak);
				}
			}
		}
		
		return peaks;
	}
}
