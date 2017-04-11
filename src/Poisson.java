import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

public class Poisson {
	
	private final double eps=1., dx = 1.;
	private int N=100;
	private int[] iPlus = new int[N+1], iMinus = new int[N+1];
	private double[][][] phi = new double[N+1][N+1][N+1], 
			phiNew = new double[N+1][N+1][N+1];
	private double[][][] rho = new double[N+1][N+1][N+1];
	private double[] E_ijk = new double[3];
	
	public Poisson(int N, double phi0) {
		this.N = N;
		this.init();
	}
	
	private void setPlusMinus(){
		iPlus = new int[N+1];
		iMinus = new int[N+1];
		for(int i=1; i<N; i++){
			iPlus[i]=i+1;
			iMinus[i]=i-1;
		}
	}
	
	private void init(){
		
		this.setPlusMinus();
		phi = new double[N+1][N+1][N+1];
		rho = new double[N+1][N+1][N+1];
		
		for(int i=1; i<N; i++)
			for(int j=1; j<N; j++)
				for(int k=1; k<N; k++){
					if (i==N/2 && j==N/2 && k==N/2) rho[i][j][k] = 1.0;
					else rho[i][j][k] = 0.;
				}
		
		this.setNew(1.);
		
	}
	
	private void getE_ijk(int i, int j, int k){
		E_ijk[0] = -(phi[iPlus[i]][j][k]-phi[iMinus[i]][j][k])/(2.*dx); 
		E_ijk[1] = -(phi[i][iPlus[j]][k]-phi[i][iMinus[j]][k])/(2.*dx); 
		E_ijk[2] = -(phi[i][j][iPlus[k]]-phi[i][j][iMinus[k]])/(2.*dx);
	}
	
	private void setNewJac(){
		
		phiNew = new double[N+1][N+1][N+1];
		
		for(int i=1; i<N; i++)
			for(int j=1; j<N; j++)
				for(int k=1; k<N; k++){
					phiNew[i][j][k] = (phi[iPlus[i]][j][k]+phi[iMinus[i]][j][k]+
							phi[i][iPlus[j]][k]+phi[i][iMinus[j]][k]+
							phi[i][j][iPlus[k]]+phi[i][j][iMinus[k]]+
							rho[i][j][k]/eps)/6.0;
				}
		
		phi = phiNew;
	}
	
	private void setNew(double omega){
		
		//phiNew = new double[N+1][N+1][N+1];
		double phi_ijk_n =0.;
		
		for(int i=1; i<N; i++)
			for(int j=1; j<N; j++)
				for(int k=1; k<N; k++){
					phi_ijk_n = phi[i][j][k];
					phi[i][j][k] = (phi[iPlus[i]][j][k]+phi[iMinus[i]][j][k]+
							phi[i][iPlus[j]][k]+phi[i][iMinus[j]][k]+
							phi[i][j][iPlus[k]]+phi[i][j][iMinus[k]]+
							rho[i][j][k]/eps)/6.0;
					phi[i][j][k] = (1.-omega)*phi_ijk_n + omega * phi[i][j][k];
				}
		
		//phi = phiNew;
	}
	
	public int evolve(String phiFile, double omega, boolean write) 
			throws FileNotFoundException, UnsupportedEncodingException{
		
		double phi_avg_temp = 0.;
		for(int i=0; i<phi.length; i++)
			for(int j=0; j<phi.length; j++)
				for(int k=0; k<phi.length; k++){
				phi_avg_temp += phi[i][j][k] / 
						(double)(N*N*N);
			}
		
		int count = 0;
		while(true){
			
			this.setNew(omega);
			count++;
			if(count%100==0 && count > 1000){
				double phi_avg = 0.;
				for(int i=1; i<N; i++)
					for (int j=1; j<N; j++)
						for(int k=1; k<N; k++)
							phi_avg += phi[i][j][k] / 
								(double)(N*N*N);
				if(Math.abs((phi_avg_temp - phi_avg)/phi_avg) < 1e-3)
					break;
				else
					phi_avg_temp = phi_avg;
			}
		}
		
		if(write){
			PrintWriter writer1 = new PrintWriter(phiFile, "UTF-8");
			for(int i=1; i<N; i++){
				for(int j=1; j<N; j++)
					for(int k=1; k<N; k++){
						this.getE_ijk(i, j, k);
						writer1.println(i + " " + j + " " + k +  " " + 
									rho[i][j][k] + " " + phi[i][j][k] + " " +
									E_ijk[0] + " " + E_ijk[1] + " " + E_ijk[2]);
					}
				writer1.println();
			}
			writer1.close();
		}
		return count;
		
	}
	
	public void SOR(String SORfile) 
			throws FileNotFoundException, UnsupportedEncodingException{
		
		PrintWriter writer1 = new PrintWriter(SORfile, "UTF-8");
		int Nomega = 50;
		double omega, d_omega = 1./(double)(Nomega);
		int nt = 0;
		for(int w=0; w<Nomega; w++){
			omega = 1. + d_omega * w;
			nt = this.evolve("", omega, false);
			writer1.println(omega + " " + nt);
			this.init();
		}
		writer1.close();
		
	}
	
	public static void main(String[] args) 
			throws FileNotFoundException, UnsupportedEncodingException {
		System.out.println("HOLA");
		Poisson poiss = new Poisson(100, 0.0);
		//poiss.evolve("Poisson_10-3.dat", 1., true);
		poiss.SOR("poissonSOR.dat");
	}
	
}
