package glycanDB;

public class Peak implements Comparable<Peak>{
	public double mz;
	public double intensity;
	public boolean annotated;

	public Peak() {
		mz = -1;
		intensity = -1;
		annotated = false;
	}
	public Peak(double mass, double intensity) {
		mz = mass;
		this.intensity = intensity;
		annotated = false;
	}

	@Override
	public int compareTo(Peak o) {
		if( o.intensity > this.intensity ) return 1;
		else if( o.intensity < this.intensity ) return -1;
		return 0;
	}
}