package glycanDB;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class GlycoSearch {
	public GlycoSearch() {
		mgf_title_offset_map = new HashMap<>();
	}
	private String mgf_path = null;
	private Map<String, Long> mgf_title_offset_map = null;
	private double TOLERANCE = 0.025;

	//	질량 상수값
	public static final double PROTON = 1.007276;
	public static final double H_2_O = 1.007825035*2 + 15.99491463;
	public static final double N_H_3 = 1.007825035*3 + 14.003074;
	
	//	아미노산의 양 끝의 (아미노산 이외의 부분의) 질량
	public static double NTermOff = PROTON;
	public static double CTermOff = H_2_O+PROTON;
	
	public static final double aminoacid_masss[]={
			71.03711,	0,			103.00919, 115.02694,	129.04259, 
			147.06841,	57.02146,	137.05891, 113.08406,	0, 
			128.09496,	113.08406,	131.04049, 114.04293,	0, 
			97.05276,	128.05858,	156.10111, 87.03203,	101.04768,
			0,			99.06841,	186.07931, 0,			163.06333,	0};
	
	
	public void readGPSM(String path) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(path));
			String line = br.readLine();	//	 header skip
			while( (line=br.readLine()) != null ) {
				String[] tok = line.split("\\t");
				String title = tok[0].trim();
				String peptide = tok[2];
				String modpeptide = tok[3];
				int cs = Integer.parseInt(tok[7]);
				double rt = Double.parseDouble(tok[8]);
				double observedM = Double.parseDouble(tok[9]);
				double cal_pep_mass = Double.parseDouble(tok[13]);
				String tmp_mod = tok[25];
				String[] glyco_id = tok[26].split("%");
				String composition = glyco_id[0].trim();
				double glycoM= Double.parseDouble(glyco_id[1].trim());

				if( mgf_title_offset_map.get(title) == null )
					continue;

				System.out.println("pep: " + peptide);
				System.out.println("mod pep: " + modpeptide);
				System.out.println("cs: " + cs);
				System.out.println("RT: " + rt);
				System.out.println("Observed Mass: " + observedM);
				System.out.println("calulated peptide mass: " + cal_pep_mass);
				System.out.println("glyco code: " + composition +"\t"+glycoM);
				System.out.println("mod :" + tmp_mod);				
				System.out.println("mod :" + tmp_mod.isEmpty());
				
				List<Peak> plist = getPeakList(mgf_path, mgf_title_offset_map.get(title));
				Set<Integer> ided_peaks = new HashSet<>();
				
//				fragment ion 계산(후보, 이론적으로 계산된 fragment ion)
				int pLen = peptide.length();
				double [] theo_b_ions = new double[pLen];
				double[] theo_y_ions = new double[pLen];
				double[] modified_position = new double[pLen];
				for(int k=0; k<pLen; k++){
					modified_position[k] = 0;
				}
				if( !tmp_mod.isEmpty() ) { 
					tmp_mod = (tmp_mod.charAt(0)=='"'?tmp_mod.substring(1):tmp_mod);
					tmp_mod = (tmp_mod.charAt(tmp_mod.length()-1)=='"'?tmp_mod.substring(0,tmp_mod.length()-1):tmp_mod);
					String[] mods = tmp_mod.split(",");
					for( int k = 0; k < mods.length; k++ ) {
						String mod = mods[k].trim();
						System.out.println(mod);
						int idxAA = mod.indexOf('(')-1;
						int pos = Integer.parseInt(mod.substring(0, idxAA));
//						char aa = mod.charAt(idxAA);
						double delta = Double.parseDouble(mod.substring(idxAA+2,mod.length()-1));
//						System.out.println("pos: "+pos + ", AA: "+aa + ", delta: " + delta);
						modified_position[pos-1] += delta;
					}
				}
				
				int i = 0;
				double tarMass = NTermOff;
				
				for( i = 0; i < pLen-1; i++){
					tarMass += aminoacid_masss[peptide.charAt(i)-'A'];
					tarMass += modified_position[i];
					theo_b_ions[i] = tarMass;
				}
				tarMass = CTermOff;
				
				for( i = pLen-1; i > 0; i-- ){	
					tarMass += aminoacid_masss[peptide.charAt(i)-'A'];
					tarMass += modified_position[i];
					theo_y_ions[i] = tarMass;
				}
				int matIndex = -1;
				for( i = theo_y_ions.length-1; i > 0; i-- ){
					matIndex = getMatchedPeak(plist, theo_y_ions[i]);
					
					if( matIndex > -1 ){
//						System.out.println("y"+ (peptide.length() - i)+"\t" + theo_y_ions[i]);
						ided_peaks.add(matIndex);
					}
						
				}

				
				
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		GlycoSearch gs = new GlycoSearch();
		
		try{
			gs.readMGF("C:\\Users\\user_lee\\Desktop\\230717\\mgf_rep1_frac1.mgf");
		}catch (IOException e) {
			e.printStackTrace();
			System.err.println(-1);
		}
		gs.readGPSM("C:\\Users\\user_lee\\Desktop\\230717\\gpsm_rep1_frac1.txt");
		
	}
	
	private void readMGF(String fpath) throws IOException{
		mgf_path = fpath;
		
		String line;
		BufferedReader br = new BufferedReader(new FileReader(fpath));
		long offset = 0;	//	charactors count
		int index = 1;
		while( (line=br.readLine()) != null ){
			offset += line.length() +2;
			if( line.startsWith("BEGIN IONS") ){
				line = br.readLine();	//	title
				String title = line.split("=")[1].trim();
				if(title.endsWith(".dta"))
					title = title.substring(0,title.length()-6);
				offset += line.length() +2;
				line = br.readLine();
				while( !line.startsWith("END IONS") ){
					String[] spliter = line.split("\\s+");
					if( spliter.length == 2 && Character.isDigit(line.charAt(0)) ){
						try{
							System.out.println(title +"\t"+offset);
							mgf_title_offset_map.put(title,offset);
							break;
						}catch(NumberFormatException e){

						}
					}

					offset += line.length() +2;
					line = br.readLine();
				}
				offset += line.length() +2;

				if( index > 11200 )
					break;
				index++;
			}
		}
		br.close();
		System.out.println(mgf_title_offset_map.size()+" scans were processed!");
		System.out.println();
	}
	
	private List<Peak> getPeakList(String file, long offset) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		List<Peak> peaklist = new ArrayList<>();
		br.skip(offset);
		String line = br.readLine();
//		System.out.println(line);
		double max_inten = -1;

		while( !line.startsWith("END IONS") ){

			String[] spliter = line.split("\\s+");
			if( spliter.length == 2 && Character.isDigit(line.charAt(0)) ){
				try{
					Peak p = new Peak();
					p.mz = Double.parseDouble(spliter[0]);
					p.intensity = Double.parseDouble(spliter[1]);
					max_inten = (max_inten<p.intensity?p.intensity:max_inten);
					peaklist.add(p);
//					System.out.println(line);
				}catch(NumberFormatException e){

				}
			}

			offset += line.length() +2;
			line = br.readLine();
		}

		br.close();
		return peaklist;
	}
	
	public int getMatchedPeak(List<Peak> peakMassList, double tarPeak){
		int matIndex = -1;
		int peaksCount = peakMassList.size();
		//	target m/z에서 허용오차값을 뺀 값의 오른쪽에 존재하는 가장 작은 값을 찾는다.
		int index = biSearchPeak(peakMassList, tarPeak - TOLERANCE, peaksCount);
		if( peakMassList.get(index).mz < tarPeak - TOLERANCE ){
			return matIndex;
		}
		
		//	target m/z값에 허용 오차를 더한 값보다 작은 범위 내에서 peak를 찾아 그 index를 반환한다.
		if( peakMassList.get(index).mz <= tarPeak + TOLERANCE ){
			matIndex = index;
		}
		return matIndex;
	}
	
	public int biSearchPeak(List<Peak> pairs, double left, int count) {
		int index;
		if( left <= pairs.get(0).mz )
			index = 0;
		else if( left > pairs.get(count-1).mz ){
			index = count-1;
		}
		else{
			int M, L =0, R = count-1;
			while( R - L > 1 ){
				M = (L + R) / 2;
				
				if( left <= pairs.get(M).mz )
					R = M;
				else
					L = M;
			}
			index = R;
		}
		return index;
	}
	
}
