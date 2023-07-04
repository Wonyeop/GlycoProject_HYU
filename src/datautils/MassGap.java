package datautils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class MassGap {
	public MassGap() {
		mass_gap = new HashMap<>();
	}
	
	public static Map<Double, Integer> mass_gap;
	
	public void print(String path) throws IOException {
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(new File(path))));
		
		out.println("Gap\tCount");
		Iterator<Double> iter = mass_gap.keySet().iterator();
		while(iter.hasNext()) {
			double gap = iter.next();
			int cnt = mass_gap.get(gap);
			out.println(gap + "\t" + cnt);
		}
		out.close();
	}
	
	public void process(String file_dir) {
		try {
			File fpath = new File(file_dir);
			if (fpath.isDirectory()) {
				File[] list_input_mgf_files = fpath.listFiles();
				for (int i = 0; i < list_input_mgf_files.length; i++) {
					System.out.println(list_input_mgf_files[i]);
					String line;
					BufferedReader br = new BufferedReader(new FileReader(list_input_mgf_files[i]));
//					int index = 1;
					while( (line=br.readLine()) != null ){
						if( line.startsWith("BEGIN IONS") ){
							line = br.readLine();
							String title= null;
							Queue_PeakList qp = null;
							String[]spliter = line.split("=");
							while( spliter.length > 1 ){
//								System.out.println("1\t"+line);
								String att = spliter[0];
								String val = spliter[1];
								if( att.compareTo("TITLE") == 0 ) {
									title = val;
								}
								line = br.readLine();	
								spliter = line.split("=");
							}
							
							while( !line.startsWith("END IONS") ){
//								System.out.println("2\t"+line);
								spliter = line.split("\\s+");
								if( spliter.length == 2 && Character.isDigit(line.charAt(0)) ){
									if( qp == null ) {
										qp = new Queue_PeakList(Double.parseDouble(spliter[0]), Double.parseDouble(spliter[1]));
									}else {
										double maxgap = qp.enQueue(new Peak(Double.parseDouble(spliter[0]), Double.parseDouble(spliter[1])));
										if( maxgap > 400 && qp.cnt > 1 ) {
											List<Double> gaps = qp.dequeue();
											for( int j = 0; j < gaps.size(); j++ ) {
												double gap = ((double)Math.round(10*gaps.get(j)))/10; 
												if( mass_gap.containsKey(gap) )
													mass_gap.put(gap, mass_gap.get(gap)+1);
												else
													mass_gap.put(gap, 1);
											}
										}
									}
								}
								line = br.readLine();
							}
							while (qp.cnt > 0) {
								List<Double> gaps = qp.dequeue();
								for (int j = 0; j < gaps.size(); j++) {
									double gap = ((double)Math.round(10*gaps.get(j)))/10; 
									if( gap > 400 ) continue;
									if (mass_gap.containsKey(gap))
										mass_gap.put(gap, mass_gap.get(gap) + 1);
									else
										mass_gap.put(gap, 1);
								}
							}
							
//							index++;
//							break;
						}
						
//						if( index > 200 )
//							break;
						
					}
					br.close();
				}
			} else {
				
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	
	// O(N^4) 
	public void process_n4(String file_dir) {
		try {
			File fpath = new File(file_dir);
			if (fpath.isDirectory()) {
				File[] list_input_mgf_files = fpath.listFiles();
				for (int i = 0; i < list_input_mgf_files.length; i++) {
					System.out.println(list_input_mgf_files[i]);
					List<Long> offsets = readMGF(list_input_mgf_files[i]);
					for( int j = 0; j < offsets.size(); j++ ) {
//						System.out.print("get peaklist");
						List<Peak> curPeaklist = getPeakList(list_input_mgf_files[i],offsets.get(j));
//						System.out.println("\t" + curPeaklist.size());
//						System.out.println("process peak list!");
						for( int k = 0; k < curPeaklist.size()-1; k++ ) {
							for( int l = k+1; l < curPeaklist.size() ; l++) {
								double m_gap =  (double)Math.round((curPeaklist.get(l).mz - curPeaklist.get(k).mz)*40)/40;
								if( m_gap > 366.140 )
									break;
								
								if( mass_gap.containsKey(m_gap) )
									mass_gap.put(m_gap, mass_gap.get(m_gap)+1);
								else
									mass_gap.put(m_gap, 1);
							}
						}
						
					}
					System.out.println("Write file!");
					print("res/" + list_input_mgf_files[i].getName().substring(0, list_input_mgf_files[i].getName().length()-4) + "_mass_gap_table.txt");
					
				}
			} else {
				List<Long> offsets = readMGF(fpath);
				for( int j = 0; j < offsets.size(); j++ ) {
					List<Peak> curPeaklist = getPeakList(fpath,offsets.get(j));
					for( int k = 0; k < curPeaklist.size()-1; k++ ) {
						for( int l = k+1; l < curPeaklist.size() ; l++) {
							double m_gap =  Math.round((curPeaklist.get(l).mz - curPeaklist.get(k).mz)*40)/40;
							if( mass_gap.containsKey(m_gap) )
								mass_gap.put(m_gap, mass_gap.get(m_gap)+1);
							else
								mass_gap.put(m_gap, 1);
						}
					}
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private List<Peak> getPeakList(File file, long offset) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		List<Peak> peaklist = new ArrayList<>();
		br.skip(offset);
		String line = br.readLine();
		
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
		List<Peak> tmp_peaks = new ArrayList<>();

//		total intensity / 100 이상인 peak		
//		for(int i = 0; i < peaklist.size(); i++) {
//			Peak p = peaklist.get(i);
//			if( p.intensity > max_inten/100 )
//				tmp_peaks.add(p);
//		}

//		intensity top 200 peak
		tmp_peaks.addAll(peaklist);
		Set<Double> top_200 = new HashSet<>();
		Collections.sort(tmp_peaks);
		for( int i = 0; i < tmp_peaks.size() ; i ++ ) {
			top_200.add(tmp_peaks.get(i).mz);
			if (i>=200) break;
		}
		tmp_peaks = new ArrayList<>();
		for (int i = 0; i < peaklist.size(); i++) {
			Peak p = peaklist.get(i);
			if (top_200.contains(p.mz))
				tmp_peaks.add(p);
		}
		
		
		br.close();
		return tmp_peaks;
	}

	private List<Long> readMGF(File fpath) throws IOException{
		String line;
		List<Long> offsets = new ArrayList<>();
		BufferedReader br = new BufferedReader(new FileReader(fpath));
		long offset = 0;	//	charactors count
//		System.out.println("Count offsets");
		int index = 1;
		while( (line=br.readLine()) != null ){
			offset += line.length() +2;
			if( line.startsWith("BEGIN IONS") ){
				line = br.readLine();
				offset += line.length() +2;
				line = br.readLine();
//				System.out.println(offset);
				
				while( !line.startsWith("END IONS") ){
					
					String[] spliter = line.split("\\s+");
					if( spliter.length == 2 && Character.isDigit(line.charAt(0)) ){
						try{
							Double.parseDouble(spliter[0]);
							offsets.add(offset);
							break;
						}catch(NumberFormatException e){
							
						}
					}
					
					offset += line.length() +2;
					line = br.readLine();
				}
//				break;
				offset += line.length() +2;
				
				index++;
			}
			
//			if( index > 200 )
//				break;
			
		}
		br.close();
		return offsets;
		
//		System.out.println("test\n\n");
//		
//		br = new BufferedReader(new FileReader(fpath));
//		br.skip(158);
//		System.out.println(br.readLine());
//		System.out.println(br.readLine());
//		br.close();
//		System.out.println("test\n\n");
//
//		br = new BufferedReader(new FileReader(fpath));
//		br.skip(5911);
//		System.out.println(br.readLine());
//		System.out.println(br.readLine());
//		br.close();
//		System.out.println("test\n\n");
//
//		br = new BufferedReader(new FileReader(fpath));
//		br.skip(11663);
//		System.out.println(br.readLine());
//		System.out.println(br.readLine());
//		br.close();
//		System.out.println("test\n\n");
//
//		br = new BufferedReader(new FileReader(fpath));
//		br.skip(17414);
//		System.out.println(br.readLine());
//		System.out.println(br.readLine());
//		br.close();
	}
	
	public static void main(String[] args) {
		MassGap mg = new MassGap();
		mg.process("D:/Glycans/PXD011533/mgf/HCDpdETXXDFT");
		
		try {
			mg.print("res/mass_gaps.txt");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
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
	
	public class Queue_PeakList {
		public Queue_PeakList() {
			cnt = 0;
			pivot = null;
			queue = new ArrayList<>();
		}
		public Queue_PeakList(double mass, double intensity) {
			pivot = new Peak(mass, intensity);
			queue = new ArrayList<>();
			enQueue(pivot);
		}
		
		
		public double enQueue(Peak p) {
			queue.add(p);
			cnt++;
			double max_gap;
			if( cnt > 1 ) {
				max_gap = p.mz - pivot.mz;
			}else {
				max_gap = p.mz;
			}
			return max_gap;
		}
		
		public List<Double> dequeue() {
			if( cnt <= 0 )	return null;
			List<Double> gaps = new ArrayList<>();
			gaps.add(pivot.mz);
			for(int i = 1; i < queue.size(); i++){
				gaps.add(queue.get(i).mz - pivot.mz);
			}
			queue.remove(0);
			cnt--;
			if( cnt > 0 ) { 
				pivot = queue.get(0);
			}
			return gaps;
		}
		
		public Peak pivot;
		public List<Peak> queue;
		public int cnt;
		
		
	}
}

