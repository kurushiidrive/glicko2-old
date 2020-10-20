import java.util.Arrays;

public class glicko2
{
	interface Compute
	{
		double func (double x);
	}
	
	public static double f (double x, Compute c)
	{
		return c.func(x);
	}
	
	public static void main (String args[])
	{
		double SYS_CONST = 0.6;																			// system constant
		double volatility = 0.06;																	// volatility measure
		double rating = 1500;																			// rating (generally from 1000-6000)
		double rd = 350;																				// rating deviation (generally from 0-350)
		int opp_rating[] = 																				// opponents' ratings
		{1391, 1391, 1665, 1665, 1665, 1665, 1237, 1237, 1674, 1674};
		int opp_rd[] = 																					// opponents' rating deviations
		{132, 132, 121, 121, 121, 121, 142, 142, 133, 133};
		int scores[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
			// win/lose
//		int wins = 0, losses = 0;
//		int periods = 7; // ???
		
		// DON'T MODIFY ANYTHING BELOW IF YOU ARE USING THIS PROGRAM;
		// IT CONTAINS THE CORE CODE FOR THE glicko-2 ALGORITHM
		
		int num_matches = scores.length;																// number of matches played
//		double win_rate = (double)wins / (wins + losses) * 100;											// win rate
		
		if (num_matches != opp_rd.length)
		{
			System.out.println("RD array invalid length");
			return;	
		}
		if (num_matches != opp_rating.length)
		{
			System.out.println("Rating array invalid length");
			return;
		}

		double volatility_p;
		double rating_p;
		double rd_p;
		double mu_p, phi_p;

		double mu = (rating-1500) / 173.7178, phi = rd / 173.7178;
		double nu, delta;
		double mu_opp[] = new double[num_matches];
		double phi_opp[] = new double[num_matches];
		for (int i = 0; i < mu_opp.length; i++)
			mu_opp[i] = (opp_rating[i]-1500) / 173.7178;
		for (int i = 0; i < mu_opp.length; i++)
			phi_opp[i] = opp_rd[i] / 173.7178;
	
		
		// set nu
		double sum1 = 0;
		for (int j = 0; j < num_matches; j++)
			sum1 += (Math.pow(g(phi_opp[j]),2) * E(mu, mu_opp[j], phi_opp[j]) * (1-E(mu, mu_opp[j], phi_opp[j])));
		nu = 1 / sum1;
		
		// set delta
		double sum2 = 0;
		for (int j = 0; j < num_matches; j++)
			sum2 += (g(phi_opp[j]) * (scores[j] - E(mu, mu_opp[j], phi_opp[j])));
		delta = nu * sum2;
		
		// get new volatility
		double a = Math.log(Math.pow(volatility, 2));
		Compute comp = (x) -> {
			return ((Math.exp(x)*(Math.pow(delta,2)-Math.pow(phi,2)-nu-Math.exp(x))) / (2*Math.pow(Math.pow(phi,2)+nu+Math.exp(x),2)))-((x-a) / Math.pow(SYS_CONST,2));
		};
		double epsilon = 0.000001;
		
		double A = a;
		double B;
		if (Math.pow(delta, 2) > Math.pow(phi, 2) + nu)
			B = Math.log(Math.pow(delta, 2)-Math.pow(phi, 2)-nu);
		else
		{
			int k = 1;
			while (f(a-k*SYS_CONST, comp) < 0)
				++k;
			B = a - k*SYS_CONST;
		}
		
		double fa = f(A, comp); double fb = f(B, comp);
		while (Math.abs(B-A) > epsilon)
		{
			double C = A + (A-B)*fa/(fb-fa);
			double fc = f(C, comp);
			if (fc*fb < 0)
			{ A = B; fa = fb; }
			else
				fa /= 2;
			B = C; fb = fc;
			if (Math.abs(B-A) <= epsilon)
				break;
		}
		
		volatility_p = Math.exp(A/2);
		
		// calibrate the rating and the rating deviation
		double phi_star = Math.sqrt(Math.pow(phi,2) + Math.pow(volatility_p,2));
		phi_p = 1 / Math.sqrt((1/Math.pow(phi_star,2)) + (1/nu));
		
		double sum3 = 0;
		for (int j = 0; j < num_matches; j++)
			sum3 += g(phi_opp[j]) * (scores[j] - E(mu, mu_opp[j], phi_opp[j]));
		mu_p = mu + Math.pow(phi_p,2)*sum3;
		
		rating_p = 173.7178 * mu_p + 1500;
		rd_p = 173.7178 * phi_p;
		
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
		
		System.out.println ("\nGlicko-2 Rating Calculator\n");
		System.out.printf ("Player Rating Before:\t%.0f\n", rating);
		System.out.printf ("Player Rating Deviation Before:\t%.0f\n", rd);
		System.out.printf ("Opponents\' Ratings:\t%s\n", Arrays.toString(opp_rating));
		System.out.printf ("Opponents\' Rating Deviations:\t%s\n", Arrays.toString(opp_rd));
		System.out.printf ("Scores against Opponents (1 win, 0 loss):\t%s\n", Arrays.toString(scores));
		System.out.printf ("Player Volatility Before:\t%.6f\n", volatility);
		System.out.printf ("System Constant:\t%.1f\n", SYS_CONST);
		System.out.printf ("Number of Matches:\t%d\n\n\n", num_matches);
		
		System.out.printf ("Player Rating Calibrated:\t%.0f\n", rating_p);
		System.out.printf ("Player Rating Deviation Calibrated:\t%.0f\n", rd_p);
		System.out.printf ("Player Volatility Calibrated:\t%.6f\n\n", volatility_p);
//		System.out.printf ("Player Record Updated:\t%d-%d\n", wins, losses);
//		System.out.printf ("Player Win Rate Updated:\t%.1f%%\n\n", win_rate);
		
		// DEBUGGING
		System.out.printf ("Nu:\t%.10f\n", nu);
		System.out.printf ("Delta:\t%.10f\n\n", delta);
	}
	
	private static double g (double p)
	{
		return 1/(Math.sqrt(1+(3*(Math.pow(p,2))/(Math.pow(Math.PI,2)))));
	}
	
	private static double E (double m, double mj, double pj)
	{
		return 1/(1+Math.exp(-1*g(pj)*(m-mj)));
	}
}
