package datautils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GapAnnotate {
	public GapAnnotate() {
		annotations = new HashMap<>();
		buildAnnotationMap();
	}
	
	public static final double PROTON = 1.007276;
	public static final double H_2_O = 1.007825035*2 + 15.99491463;
	public static final double N_H_3 = 1.007825035*3 + 14.003074;
	
	public static final double glycans[] = {
						// 	Agg.-NLoss-NLoss:	Residue Name				comp.		mono		avg.
			132.0423, 	// 	-NLoss-NLoss: 		Pentose						C5O4H8 		132.0423	132.12
			146.0579, 	// 	F-NLoss-NLoss: 		Fucose(Deoxyhexose)			C6O4H10 	146.0579 	146.14
			162.0528, 	// 	H-NLoss-NLoss: 		Hexose(Mannose,Galactose) 	C6O5H10 	162.0528 	162.14
			161.0688,	//	-NLoss-NLoss:		Hexosamine					C6O4NH11	161.0688	161.16
			176.0321,	//	-NLoss-NLoss:		Hexuronic acid				C6O6H8		176.0321	176.13
			198.0140,	//	Na1-NLoss-NLoss:	Na salt 					C6O6H7Na	198.0140	198.11
			192.0634,	//	-NLoss-NLoss:		Heptose						C7O6H12		192.0634	192.17
			203.0794, 	// 	N-NLoss-NLoss: 		N-Acetylaminohexose 		C8O5NH13 	203.0794 	203.19	(HexNAc)
			220.0583, 	// 	-NLoss-NLoss:		KDO							C8O7H12		220.0583	220.18
			291.0954,	//	S-NLoss-NLoss:		N-Acetylneuraminic acid		C11O8NH17	291.0954	291.26	Sialic acid(Neu5Ac)
			313.0773,	//	Na2-NLoss-NLoss:	Na salt 					C11O8NH16Na	313.0773	313.24
			307.0903,	//	-NLoss-NLoss:		N-Glycoylneuraminic acid	C11O9NH17	307.0903	307.26	(Neu5Gc)
			329.0734,	//	Na3-NLoss-NLoss:	Na salt 					C11O9NH16Na	329.0734	329.24
	};
	//

	
	public static final double oxoniums[] = { 
			126.055-PROTON,	//	C6H7NO2+	[HexNAc - C2H6O3] 78.032
			138.055-PROTON,	//	C7H8NO2+	[HexNAc - CH6O3]66.032
			144.065-PROTON,	//	C6H10NO3+	[HexNAc - C2H4O2]60.022	
			168.066-PROTON,	//	C8H10NO3+	[HexNAc - 2H2O]36.021
			186.076-PROTON,	//	C5H12NO4+	[HexNAc - H2O]18.011
//			204.087,	//	C8H14NO5+	[HexNAc]
			274.092-PROTON,	//	C11H6NO7+	[Neu5Ac-H2O]
//			292.103,	//	C11H18NO8+	[Neu5Ac]
//			366.140-PROTON		//	C14H24NO10+	[HexHexNAc]
	};
	
	public static final double aminoacid_masss[]={
			71.03711,	0,			103.00919, 115.02694,	129.04259,		//ABCDE 
			147.06841,	57.02146,	137.05891, 113.08406,	0, 				//FGHIJ
			128.09496,	113.08406,	131.04049, 114.04293,	0, 				//KLMNO
			97.05276,	128.05858,	156.10111, 87.03203,	101.04768,		//PQRST
			0,			99.06841,	186.07931, 0,			163.06333,	0};	//UVWXY
	
	public Map<Double, String> annotations = null;
	
	public List<Double> masses;
	
	private void buildAnnotationMap() {
		for( int m = 1; m < 5; m++ ) {
			annotations.put(((double)(Math.round(aminoacid_masss['A'-'A'])*10)/(m*10)), "AA A,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['C'-'A'])*10)/(m*10)), "AA C,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['D'-'A'])*10)/(m*10)), "AA D,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['E'-'A'])*10)/(m*10)), "AA E,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['F'-'A'])*10)/(m*10)), "AA F,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['G'-'A'])*10)/(m*10)), "AA G,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['H'-'A'])*10)/(m*10)), "AA H,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['I'-'A'])*10)/(m*10)), "AA I,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['K'-'A'])*10)/(m*10)), "AA K,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['L'-'A'])*10)/(m*10)), "AA L,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['M'-'A'])*10)/(m*10)), "AA M,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['N'-'A'])*10)/(m*10)), "AA N,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['P'-'A'])*10)/(m*10)), "AA P,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['Q'-'A'])*10)/(m*10)), "AA Q,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['R'-'A'])*10)/(m*10)), "AA R,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['S'-'A'])*10)/(m*10)), "AA S,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['T'-'A'])*10)/(m*10)), "AA T,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['V'-'A'])*10)/(m*10)), "AA V,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['W'-'A'])*10)/(m*10)), "AA W,CS "+m);
			annotations.put(((double)(Math.round(aminoacid_masss['Y'-'A'])*10)/(m*10)), "AA Y,CS "+m);
			
			annotations.put(((double)(Math.round(oxoniums[0])*10)/(m*10)), "oxonium: HexNAc - C2H6O3,CS "+m);
			annotations.put(((double)(Math.round(oxoniums[1])*10)/(m*10)), "oxonium: HexNAc - CH6O3,CS "+m);
			annotations.put(((double)(Math.round(oxoniums[2])*10)/(m*10)), "oxonium: HexNAc - C2H4O2,CS "+m);
			annotations.put(((double)(Math.round(oxoniums[3])*10)/(m*10)), "oxonium: HexNAc - 2H2O,CS "+m);
			annotations.put(((double)(Math.round(oxoniums[4])*10)/(m*10)), "oxonium: HexNAc - H2O,CS "+m);
			annotations.put(((double)(Math.round(oxoniums[5])*10)/(m*10)), "oxonium: Neu5Ac-H2O,CS "+m);
			
			annotations.put(((double)(Math.round(glycans[0])*10)/(m*10)), "glycan: Pentose,CS "+m);
			annotations.put(((double)(Math.round(glycans[1])*10)/(m*10)), "glycan: Fucose,CS "+m);
			annotations.put(((double)(Math.round(glycans[2])*10)/(m*10)), "glycan: Hexose,CS "+m);
			annotations.put(((double)(Math.round(glycans[3])*10)/(m*10)), "glycan: Hexosamine,CS "+m);
			annotations.put(((double)(Math.round(glycans[4])*10)/(m*10)), "glycan: Hexuronic acid,CS "+m);
			annotations.put(((double)(Math.round(glycans[5])*10)/(m*10)), "glycan: Na salt,CS "+m);
			annotations.put(((double)(Math.round(glycans[6])*10)/(m*10)), "glycan: Heptose,CS "+m);
			annotations.put(((double)(Math.round(glycans[7])*10)/(m*10)), "glycan: HexNAc,CS "+m);
			annotations.put(((double)(Math.round(glycans[8])*10)/(m*10)), "glycan: KDO,CS "+m);
			annotations.put(((double)(Math.round(glycans[9])*10)/(m*10)), "glycan: Neu5Ac,CS "+m);
			annotations.put(((double)(Math.round(glycans[10])*10)/(m*10)), "glycan: Na salt ,CS "+m);
			annotations.put(((double)(Math.round(glycans[11])*10)/(m*10)), "glycan: Neu5Gc,CS "+m);
			annotations.put(((double)(Math.round(glycans[12])*10)/(m*10)), "glycan: Na salt ,CS "+m);

			annotations.put(((double)(Math.round(H_2_O)*10)/(m*10)), "NLoss: H2O,CS "+m);
			annotations.put(((double)(Math.round(N_H_3)*10)/(m*10)), "NLoss: NH3,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['A'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: A,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['C'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: C,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['D'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: D,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['E'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: E,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['F'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: F,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['G'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: G,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['H'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: H,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['I'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: I,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['K'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: K,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['L'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: L,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['M'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: M,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['N'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: N,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['P'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: P,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['Q'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: Q,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['R'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: R,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['S'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: S,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['T'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: T,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['V'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: V,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['W'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: W,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['Y'-'A']-H_2_O)*10)/(m*10)), "AA-NLoss: Y,CS "+m);
//			
//			annotations.put(((double)(Math.round(oxoniums[0]-H_2_O)*10)/(m*10)), "oxonium-NLoss: HexNAc - C2H6O3,CS "+m);
//			annotations.put(((double)(Math.round(oxoniums[1]-H_2_O)*10)/(m*10)), "oxonium-NLoss: HexNAc - CH6O3,CS "+m);
//			annotations.put(((double)(Math.round(oxoniums[2]-H_2_O)*10)/(m*10)), "oxonium-NLoss: HexNAc - C2H4O2,CS "+m);
//			annotations.put(((double)(Math.round(oxoniums[3]-H_2_O)*10)/(m*10)), "oxonium-NLoss: HexNAc - 2H2O,CS "+m);
//			annotations.put(((double)(Math.round(oxoniums[4]-H_2_O)*10)/(m*10)), "oxonium-NLoss: HexNAc - H2O,CS "+m);
//			annotations.put(((double)(Math.round(oxoniums[5]-H_2_O)*10)/(m*10)), "oxonium-NLoss: Neu5Ac-H2O,CS "+m);
//			
//			annotations.put(((double)(Math.round(glycans[0]-H_2_O)*10)/(m*10)), "glycan-NLoss: Pentose,CS "+m);
//			annotations.put(((double)(Math.round(glycans[1]-H_2_O)*10)/(m*10)), "glycan-NLoss: Fucose,CS "+m);
//			annotations.put(((double)(Math.round(glycans[2]-H_2_O)*10)/(m*10)), "glycan-NLoss: Hexose,CS "+m);
//			annotations.put(((double)(Math.round(glycans[3]-H_2_O)*10)/(m*10)), "glycan-NLoss: Hexosamine,CS "+m);
//			annotations.put(((double)(Math.round(glycans[4]-H_2_O)*10)/(m*10)), "glycan-NLoss: Hexuronic acid,CS "+m);
//			annotations.put(((double)(Math.round(glycans[5]-H_2_O)*10)/(m*10)), "glycan-NLoss: Na salt,CS "+m);
//			annotations.put(((double)(Math.round(glycans[6]-H_2_O)*10)/(m*10)), "glycan-NLoss: Heptose,CS "+m);
//			annotations.put(((double)(Math.round(glycans[7]-H_2_O)*10)/(m*10)), "glycan-NLoss: HexNAc,CS "+m);
//			annotations.put(((double)(Math.round(glycans[8]-H_2_O)*10)/(m*10)), "glycan-NLoss: KDO,CS "+m);
//			annotations.put(((double)(Math.round(glycans[9]-H_2_O)*10)/(m*10)), "glycan-NLoss: Neu5Ac,CS "+m);
//			annotations.put(((double)(Math.round(glycans[10]-H_2_O)*10)/(m*10)), "glycan-NLoss: Na salt ,CS "+m);
//			annotations.put(((double)(Math.round(glycans[11]-H_2_O)*10)/(m*10)), "glycan-NLoss: Neu5Gc,CS "+m);
//			annotations.put(((double)(Math.round(glycans[12]-H_2_O)*10)/(m*10)), "glycan-NLoss: Na salt ,CS "+m);
//			
//			annotations.put(((double)(Math.round(aminoacid_masss['A'-'A']-N_H_3)*10)/(m*10)), "AA A,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['C'-'A']-N_H_3)*10)/(m*10)), "AA C,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['D'-'A']-N_H_3)*10)/(m*10)), "AA D,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['E'-'A']-N_H_3)*10)/(m*10)), "AA E,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['F'-'A']-N_H_3)*10)/(m*10)), "AA F,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['G'-'A']-N_H_3)*10)/(m*10)), "AA G,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['H'-'A']-N_H_3)*10)/(m*10)), "AA H,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['I'-'A']-N_H_3)*10)/(m*10)), "AA I,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['K'-'A']-N_H_3)*10)/(m*10)), "AA K,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['L'-'A']-N_H_3)*10)/(m*10)), "AA L,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['M'-'A']-N_H_3)*10)/(m*10)), "AA M,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['N'-'A']-N_H_3)*10)/(m*10)), "AA N,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['P'-'A']-N_H_3)*10)/(m*10)), "AA P,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['Q'-'A']-N_H_3)*10)/(m*10)), "AA Q,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['R'-'A']-N_H_3)*10)/(m*10)), "AA R,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['S'-'A']-N_H_3)*10)/(m*10)), "AA S,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['T'-'A']-N_H_3)*10)/(m*10)), "AA T,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['V'-'A']-N_H_3)*10)/(m*10)), "AA V,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['W'-'A']-N_H_3)*10)/(m*10)), "AA W,CS "+m);
//			annotations.put(((double)(Math.round(aminoacid_masss['Y'-'A']-N_H_3)*10)/(m*10)), "AA Y,CS "+m);
//			
//			annotations.put(((double)(Math.round(oxoniums[0]-N_H_3)*10)/(m*10)), "oxonium-NLoss: HexNAc - C2H6O3,CS "+m);
//			annotations.put(((double)(Math.round(oxoniums[1]-N_H_3)*10)/(m*10)), "oxonium-NLoss: HexNAc - CH6O3,CS "+m);
//			annotations.put(((double)(Math.round(oxoniums[2]-N_H_3)*10)/(m*10)), "oxonium-NLoss: HexNAc - C2H4O2,CS "+m);
//			annotations.put(((double)(Math.round(oxoniums[3]-N_H_3)*10)/(m*10)), "oxonium-NLoss: HexNAc - 2H2O,CS "+m);
//			annotations.put(((double)(Math.round(oxoniums[4]-N_H_3)*10)/(m*10)), "oxonium-NLoss: HexNAc - H2O,CS "+m);
//			annotations.put(((double)(Math.round(oxoniums[5]-N_H_3)*10)/(m*10)), "oxonium-NLoss: Neu5Ac-H2O,CS "+m);
//			
//			annotations.put(((double)(Math.round(glycans[0]-N_H_3)*10)/(m*10)), "glycan-NLoss: Pentose,CS "+m);
//			annotations.put(((double)(Math.round(glycans[1]-N_H_3)*10)/(m*10)), "glycan-NLoss: Fucose,CS "+m);
//			annotations.put(((double)(Math.round(glycans[2]-N_H_3)*10)/(m*10)), "glycan-NLoss: Hexose,CS "+m);
//			annotations.put(((double)(Math.round(glycans[3]-N_H_3)*10)/(m*10)), "glycan-NLoss: Hexosamine,CS "+m);
//			annotations.put(((double)(Math.round(glycans[4]-N_H_3)*10)/(m*10)), "glycan-NLoss: Hexuronic acid,CS "+m);
//			annotations.put(((double)(Math.round(glycans[5]-N_H_3)*10)/(m*10)), "glycan-NLoss: Na salt,CS "+m);
//			annotations.put(((double)(Math.round(glycans[6]-N_H_3)*10)/(m*10)), "glycan-NLoss: Heptose,CS "+m);
//			annotations.put(((double)(Math.round(glycans[7]-N_H_3)*10)/(m*10)), "glycan-NLoss: HexNAc,CS "+m);
//			annotations.put(((double)(Math.round(glycans[8]-N_H_3)*10)/(m*10)), "glycan-NLoss: KDO,CS "+m);
//			annotations.put(((double)(Math.round(glycans[9]-N_H_3)*10)/(m*10)), "glycan-NLoss: Neu5Ac,CS "+m);
//			annotations.put(((double)(Math.round(glycans[10]-N_H_3)*10)/(m*10)), "glycan-NLoss: Na salt ,CS "+m);
//			annotations.put(((double)(Math.round(glycans[11]-N_H_3)*10)/(m*10)), "glycan-NLoss: Neu5Gc,CS "+m);
//			annotations.put(((double)(Math.round(glycans[12]-N_H_3)*10)/(m*10)), "glycan-NLoss: Na salt ,CS "+m);
			
			
		}
		annotations.put(1.0, "Isotope");
		annotations.put(2.0, "Isotope");
		annotations.put(3.0, "Isotope");
		annotations.put(0.5, "Isotope cs2");
		annotations.put(0.3, "Isotope cs3");
		
		masses = new ArrayList<>();
		masses.addAll(annotations.keySet());
		Collections.sort(masses);
		
//		Iterator<Double> iter = annotations.keySet().iterator();
//		while( iter.hasNext() ) System.out.println(iter.next());
	}
	
	public int biSearchPeak(List<Double> pairs, double left, int count) {
		int index;
		if( left <= pairs.get(0) )
			index = 0;
		else if( left > pairs.get(count-1) ){
			index = count-1;
		}
		else{
			int M, L =0, R = count-1;
			while( R - L > 1 ){
				M = (L + R) / 2;
				
				if( left <= pairs.get(M) )
					R = M;
				else
					L = M;
			}
			index = R;
		}
		return index;
	}
	private double TOLERANCE = 0.05;
	public int getMatchedPeak(List<Double> peakMassList, double tarPeak, double minimum){
//		double maxPeak = -1;//, maxDelta = TOLERANCE+1;
		int matIndex = -1;
		int peaksCount = peakMassList.size();
		//	target m/z에서 허용오차값을 뺀 값의 오른쪽에 존재하는 가장 작은 값을 찾는다.
		int index = biSearchPeak(peakMassList, tarPeak - TOLERANCE, peaksCount);
//		System.out.println(spectrum.peakMassList.size() + "  " + index + "  " + tarPeak + "\t" + spectrum.peakMassList.get(index) + "\t" + spectrum.peakMassList.get(19));
		if( peakMassList.get(index) < tarPeak - TOLERANCE ){
			return matIndex;
			//return index;
		}
		//int index = 0;
		
		
		//	target m/z값에 허용 오차를 더한 값보다 작은 범위 내에서 peak를 찾아 그 index를 반환한다.
		if( peakMassList.get(index) <= tarPeak + TOLERANCE ){
			matIndex = index;
		}
		return matIndex;
	}
	
	public void process(String path) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(path));
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(new File(path.substring(0,path.length()-4)+"_anno.txt" ))));
		String line=br.readLine();	//	header
		out.println(line + "\tAnnotation");
		while( (line=br.readLine()) != null ) {
			String[] tok = line.split("\t");
			double gap = Double.parseDouble(tok[0]);
			
			int idx = getMatchedPeak(masses, gap, 0);
			
			if( idx < 0 )	{
				out.println(line+"\tNA");
				continue;
			}
			String anno = annotations.get(masses.get(idx));
			out.println(line+"\t"+anno);
		}
		br.close();
		out.close();		
	}
	
	public static void main(String[] args) {
		GapAnnotate GA = new GapAnnotate();
		System.out.println(GA.annotations.size());
		try {
			GA.process("res/mass_gaps.txt");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
}
