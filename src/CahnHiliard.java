
public class CahnHiliard {
	
	private final double a=0.1, M=0.1, kappa=0.01, noise=0.001,
			dx = 0.1, dt = 0.01;
	private int N=100;
	private int[] iPlus = new int[N], iMinus = new int[N];
	private double[][] phi = new double[N][N], mu = new double[N][N];
	private double[][] phiNew = new double[N][N], muNew = new double[N][N];
	
	public CahnHiliard(int N, double phi0) {
		this.N = N;
		this.init(phi0);
	}
	
	private void setPlusMinus(){
		for(int i=0; i<N; i++){
			iPlus[i]=i+1;
			iMinus[i]=i-1;
		}
		iPlus[N-1]=0;
		iMinus[0]=N-1;
	}
	
	private void setMu(){
		double phi_ij = phi[0][0], d2phi_ij=0;
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++){
				phi_ij = phi[i][j];
				d2phi_ij = phi[iPlus[i]][j]+phi[iMinus[i]][j]+phi[i][iPlus[j]]+
						phi[i][iMinus[j]]-4*phi_ij;
				d2phi_ij/=(dx*dx);
				mu[i][j] = a*phi_ij+a*phi_ij*phi_ij*phi_ij-kappa*d2phi_ij;
			}
	}
	
	private void init(double phi0){
		this.setPlusMinus();
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++)
				phi[i][j] = phi0 + (Math.random()-0.5)*noise;
		this.setMu();
	}
	
	private void setNew(){
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++){
				phiNew[i][j] = phi[i][j] + M*dt*(mu[iPlus[i]][j]+mu[iMinus[i]][j]+
						mu[i][iPlus[j]]+mu[i][iMinus[j]]-4*mu[i][j])/(dx*dx);
			}
		phi = phiNew;
		this.setMu();
	}
	
}
