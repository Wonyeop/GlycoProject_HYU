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
	
	public static final double NEUTRON = 1.008665;  
	private double TOLERANCE = 0.025;

	public static Map<Double, Integer> mass_gap;
	
	public void print(String path) {
		try {
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(new File(path))));
			
			out.println("Gap\tCount");
			Iterator<Double> iter = mass_gap.keySet().iterator();
			while(iter.hasNext()) {
				double gap = iter.next();
				int cnt = mass_gap.get(gap);
				out.println(gap + "\t" + cnt);
			}
			out.close();			
		}catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void process(String file_dir) {
		try {
			File fpath = new File(file_dir);
			if (fpath.isDirectory()) {
				System.out.println("Run!");
				File[] list_input_mgf_files = fpath.listFiles();
				for (int i = 0; i < list_input_mgf_files.length; i++) {
					processFile(list_input_mgf_files[i]);
				}
			} else {
				processFile(fpath);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

	private void processFile(File mgf_files) throws FileNotFoundException, IOException {
		System.out.println(mgf_files);
		String line;
		BufferedReader br = new BufferedReader(new FileReader(mgf_files));
//					int index = 1;
		while( (line=br.readLine()) != null ){
			if( line.startsWith("BEGIN IONS") ){
				line = br.readLine();
				String title= null;
				Queue_PeakList qp = null;
				String[]spliter = line.split("=");
				
//							skipping meta info.
				while( spliter.length > 1 ){
//					System.out.println("1\t"+line);
					String att = spliter[0];
					String val = spliter[1];
					if( att.compareTo("TITLE") == 0 ) {
						title = val;
					}
					line = br.readLine();	
					spliter = line.split("=");
				}
				
				while( !line.startsWith("END IONS") ){
//					System.out.println("2\t"+line);
					spliter = line.split("\\s+");
					if( spliter.length == 2 && Character.isDigit(line.charAt(0)) ){
						Peak p = new Peak(Double.parseDouble(spliter[0]), Double.parseDouble(spliter[1]), title);
						if( qp == null ) {
							qp = new Queue_PeakList(p);
						}else {
							double maxgap = qp.enQueue(p);
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
			}
		}
		br.close();
	}

	public static void main(String[] args) {
		long beforeTime = System.currentTimeMillis(); // 코드 실행 전에 시간 받아오기

		MassGap mg = new MassGap();
		mg.process("D:/Glycans/PXD011533/mgf/HCDpdETXXDFT/");

		mg.print("res/mass_gaps.txt");
		
		long afterTime = System.currentTimeMillis(); // 코드 실행 후에 시간 받아오기
		long diffTime = afterTime - beforeTime; // 두 개의 실행 시간
		System.out.println("실행 시간(ms): " + diffTime); // 세컨드(초 단위 변환)
	}
	
	public class Peak implements Comparable<Peak>{
		public double mz;
		public double intensity;
		public int cs;
		public boolean deconv;
		public boolean annotated;
		public String title;
		
		public Peak() {
			deconv = false;
			mz = -1;
			intensity = -1;
			cs = -1;
			annotated = false;
		}
		public Peak(double mass, double intensity, String title) {
			mz = mass;
			this.intensity = intensity;
			annotated = false;
			deconv = false;
			cs = 1;
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
		public Queue_PeakList(Peak p) {
			pivot = p;
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
		
		
		public List<Peak> deconv(List<Peak> process_queue){
			if( process_queue.size() < 2 )
				return process_queue;
			List<Peak> result_queue = new ArrayList<>();
			List<Peak> tmp_queue = new ArrayList<>();

			Peak target = process_queue.get(0);
			double prev_mz = target.mz;
			boolean chk_iso = false;
			int cs = 4;
			for (int i = 1; i < process_queue.size(); i++) {
				Peak tmp_p = process_queue.get(i);
				double gap = tmp_p.mz - prev_mz;
				if( chk_iso ) {
					if( gap > (NEUTRON/cs) + TOLERANCE ) {
						for( ; i < process_queue.size(); i++ )
							tmp_queue.add(process_queue.get(i));
						break;
					}
					if( (NEUTRON/cs) - TOLERANCE <= gap &&  gap <= (NEUTRON/cs) + TOLERANCE ) {
						prev_mz = target.mz;
					}
				}else {
					for (; cs > 0; cs--) {
						if( (NEUTRON/cs) - TOLERANCE <= gap &&  gap <= (NEUTRON/cs) + TOLERANCE ) {
							prev_mz = target.mz;
							target.cs = cs;
							chk_iso = true;
							break;
						}
					}
					if( !chk_iso ) tmp_queue.add(tmp_p);
				}
			}	
			result_queue.add(target);
			result_queue.addAll(deconv(tmp_queue));
			return result_queue;
		}
		
		public List<Double> dequeue() {
			List<Peak> processsed_queue = deconv(queue);
			if( cnt <= 0 )	return null;
			List<Double> gaps = new ArrayList<>();
			gaps.add(pivot.mz);
			for(int i = 1; i < processsed_queue.size(); i++){
				Peak tmpp = processsed_queue.get(i);
				if( pivot.cs == tmpp.cs )
					gaps.add(tmpp.mz- pivot.mz);
				
			}
			processsed_queue.remove(0);
			cnt--;
			if( cnt > 0 ) { 
				pivot = processsed_queue.get(0);
			}
			return gaps;
		}
		
		public Peak pivot;
		public List<Peak> queue;
		public int cnt;
		
		
	}
}

