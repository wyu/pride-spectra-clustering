package uk.ac.ebi.pride.tools.pride_spectra_clustering.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.ParamGroup;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.SpectrumImplementation;

public class ClusteringSpectrum {
	private String id;
	private Double precursorMZ;
	private Double precursorIntensity;
	private Integer precursorCharge;
	private List<Peak> peaklist;
	private ParamGroup additional;
	private Integer msLevel;
	
	public ClusteringSpectrum(Spectrum s) {
		id = s.getId();
		precursorMZ = s.getPrecursorMZ();
		precursorIntensity = s.getPrecursorIntensity();
		precursorCharge = s.getPrecursorCharge();
		additional = s.getAdditional();
		msLevel = s.getMsLevel();
		
		peaklist = new ArrayList<Peak>(s.getPeakList().size());
		
		for (Double mz : s.getPeakList().keySet())
			peaklist.add(new Peak(mz, s.getPeakList().get(mz)));
		
		Collections.sort(peaklist, PeakIntensityComparator.getInstance());
	}
	
	public ClusteringSpectrum(String id, Double precursorMZ,
			Double precursorIntensity, Integer precursorCharge,
			Map<Double, Double> peaklist, ParamGroup additional, Integer msLevel) {
		this.id = id;
		this.precursorMZ = precursorMZ;
		this.precursorIntensity = precursorIntensity;
		this.precursorCharge = precursorCharge;
		this.additional = additional;
		this.msLevel = msLevel;
		
		this.peaklist = new ArrayList<Peak>(peaklist.size());
		
		for (Double mz : peaklist.keySet())
			this.peaklist.add(new Peak(mz, peaklist.get(mz)));
		
		Collections.sort(this.peaklist, PeakIntensityComparator.getInstance());
	}



	public String getId() {
		return id;
	}
	public Double getPrecursorMZ() {
		return precursorMZ;
	}
	public Double getPrecursorIntensity() {
		return precursorIntensity;
	}
	public Integer getPrecursorCharge() {
		return precursorCharge;
	}
	public List<Peak> getPeaklist() {
		return peaklist;
	}
	public ParamGroup getAdditional() {
		return additional;
	}
	
	public Spectrum toSpectrum() {
		Map<Double, Double> mapPeaks = new HashMap<Double, Double>(peaklist.size());
		for (Peak p : peaklist)
			mapPeaks.put(p.getMz(), p.getIntensity());
		return new SpectrumImplementation(id, precursorCharge, precursorMZ, precursorIntensity, mapPeaks, msLevel);
	}
}
