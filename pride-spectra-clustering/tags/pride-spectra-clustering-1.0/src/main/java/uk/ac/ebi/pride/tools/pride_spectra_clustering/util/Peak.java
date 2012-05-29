package uk.ac.ebi.pride.tools.pride_spectra_clustering.util;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class Peak {
	private final double mz;
	private final double intensity;
	/**
	 * Number of times the peak was observed
	 */
	private final int count;
	
	public Peak(Peak peak) {
		this.mz = peak.mz;
		this.intensity = peak.intensity;
		this.count = peak.count;
	}
	
	public Peak(double mz, double intensity) {
		this.mz = mz;
		this.intensity = intensity;
		this.count = 1;
	}
	public Peak(double mz, double intensity, int count) {
		this.mz = mz;
		this.intensity = intensity;
		this.count = count;
	}
	public double getMz() {
		return mz;
	}
	public double getIntensity() {
		return intensity;
	}
	public int getCount() {
		return count;
	}
	@Override
	public String toString() {
		return "{m/z = " + mz + ", intensity = " + intensity + ", count = " + count + "}";
	}
	
	/**
	 * Converts a peak list stored as a Map
	 * with the m/z value as key and the
	 * intensities as values to a List of
	 * PeakS.
	 * @param peakList
	 * @return
	 */
	public static List<Peak> peakListFromMap(Map<Double, Double> peakList) {
		if (peakList == null)
			return null;
		
		List<Peak> peaks = new ArrayList<Peak>(peakList.size());
		
		for (Double mz : peakList.keySet())
			peaks.add(new Peak(mz, peakList.get(mz)));
		
		return peaks;
	}
}
