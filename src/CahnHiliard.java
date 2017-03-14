import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

public class CahnHiliard {
	
	private final double a=0.1, M=0.1, kappa=0.01, noise=0.01,
			dx = 10., dt = 0.01;
	private int N=100;
	private int[] iPlus = new int[N], iMinus = new int[N];
	private double[][] phi = new double[N][N], mu = new double[N][N], f = new double[N][N];
	private double[][] phiNew = new double[N][N], muNew = new double[N][N], fNew = new double[N][N] ;
	
	public CahnHiliard(int N, double phi0) {
		this.N = N;
		this.init(phi0);
	}
	
	private void setPlusMinus(){
		iPlus = new int[N];
		iMinus = new int[N];
		for(int i=0; i<N; i++){
			iPlus[i]=i+1;
			iMinus[i]=i-1;
		}
		iPlus[N-1]=0;
		iMinus[0]=N-1;
	}
	
	private void setMu(){
		mu = new double[N][N];
		muNew = new double[N][N];
		double phi_ij = 0., d2phi_ij=0.;
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++){
				phi_ij = phi[i][j];
				d2phi_ij = phi[iPlus[i]][j]+phi[iMinus[i]][j]+phi[i][iPlus[j]]+
						phi[i][iMinus[j]]-4*phi_ij;
				d2phi_ij/=(dx*dx);
				mu[i][j] = -a*phi_ij+a*phi_ij*phi_ij*phi_ij-kappa*d2phi_ij;
			}
	}
	
	private void setFree(){
		f = new double[N][N];
		fNew = new double[N][N];
		double phi_ij = 0., dphi_ij=0.;
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++){
				phi_ij = phi[i][j];
				dphi_ij = phi[iPlus[i]][j]+phi[iMinus[i]][j]+phi[i][iPlus[j]]+phi[i][iMinus[j]]-4*phi_ij;
				dphi_ij/=dx;
				f[i][j] = -a*phi_ij*phi_ij/2.+a*phi_ij*phi_ij*phi_ij*phi_ij/4.+kappa*dphi_ij*dphi_ij/2.;
			}
	}
	
	private void init(double phi0){
		phi = new double[N][N];
		phiNew = new double[N][N];
		this.setPlusMinus();
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++)
				phi[i][j] = phi0 + (Math.random()-0.5)*noise;
		this.setMu();
		this.setFree();
	}
	
	private void setNew(){
		
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++){
				phiNew[i][j] = phi[i][j] + M*dt*(mu[iPlus[i]][j]+mu[iMinus[i]][j]+mu[i][iPlus[j]]+
						mu[i][iMinus[j]]-4*mu[i][j])/(dx*dx);
			}
		phi = phiNew;
		this.setMu();
		this.setFree();
	}
	
	private double totalFree(){
		double free = 0.;
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++)
				free += f[i][j];
		return free;
	}
	
	public void evolve(String outFile, String freeFile, int dataPoints, int nt) 
			throws FileNotFoundException, UnsupportedEncodingException{
		
		PrintWriter writer2 = new PrintWriter(freeFile, "UTF-8");
		for(int t=0; t<nt; t++){
			if(t%dataPoints==0 && t!=0){
				PrintWriter writer1 = new PrintWriter(outFile, "UTF-8");
				writer1.printf(N + " ");
				for(int j=0; j<N; j++) writer1.printf(j + " ");
				writer1.println();
				for(int i=0; i<N; i++){
					writer1.printf(i + " ");
					for(int j=0; j<N; j++)
						writer1.printf(phi[i][j] + " ");
					writer1.println();
				}
				writer1.close();
				writer2.println((int)(t/dataPoints) + " " + this.totalFree());
			}
			this.setNew();
		}
		writer2.close();
		
	}
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
		System.out.println("HOLA");
		CahnHiliard CH = new CahnHiliard(50, -0.5);
		int dataPoints = 1000;
		CH.evolve("CHphiDensity.dat", "CHfreeEnergyDensity.dat", dataPoints, 10000000);
	}
	
}
