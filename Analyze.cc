#include <iostream>
#include <math.h>

struct pop_values {
  // Females without preference
  double female_pTT = 1.0; // f1
  double female_ptT = 1.0; // f2
  double female_ptt = 1.0; // f3

  // Females with preference
  double female_PTT = 1.0; // f4
  double female_PtT = 1.0; // f5
  double female_Ptt = 1.0; // f6

  // Males
  double male_TT = 1.0;    // f7
  double male_tT = 1.0;    // f8
  double male_tt = 1.0;    // f9

  pop_values() { Normalize(); }
  pop_values(const pop_values &) = default;
  ~pop_values() { ; }
  pop_values & operator=(const pop_values &) = default;

  void Print(std::ostream & os = std::cout) {
    os << "female_pTT = " << female_pTT << "   ";
    os << "female_ptT = " << female_ptT << "   ";
    os << "female_ptt = " << female_ptt << "\n";

    // Females with preference
    os << "female_PTT = " << female_PTT << "   ";
    os << "female_PtT = " << female_PtT << "   ";
    os << "female_Ptt = " << female_Ptt << "\n";

    // Males
    os << "male_TT = " << male_TT << "      ";
    os << "male_tT = " << male_tT << "      ";
    os << "male_tt = " << male_tt << "\n";
    
    // Genes
    os << "P = " << (female_PTT + female_PtT + female_Ptt) << "  ";
    os << "T = " << (male_TT + male_tT / 2) << "  ";
    os << "t = " << (male_tt + male_tT / 2) << "\n";
  }

  void MinPrint(size_t ud, std::ostream & os = std::cout) {
    os << ud << ", "
       << (female_PTT + female_PtT + female_Ptt) << ", "
       << male_TT << ", "
       << male_tT << ", "
       << male_tt << std::endl;
    }

  void Normalize() {
    double f_tot = female_pTT + female_ptT + female_ptt + female_PTT + female_PtT + female_Ptt;
    female_pTT /= f_tot;
    female_ptT /= f_tot;
    female_ptt /= f_tot;
    female_PTT /= f_tot;
    female_PtT /= f_tot;
    female_Ptt /= f_tot;

    male_TT = female_pTT + female_PTT;
    male_tT = female_ptT + female_PtT;
    male_tt = female_ptt + female_Ptt;
  }

  void Run(
		  const size_t updates = 1,         // How long should this run for?
		  const double female_gain = 0.05,  // s_f
		  const double male_harm = 0.5,     // s_m
		  const double pref_level = 7.5,    // alpha
		  const double dominance = 0.5,     // h_T  (and h_f?)
      const double t_mut_rate = 0.001,  // u ; Mutation rate t <-> T
      const double p_mut_rate = 0.0     // not in paper ; Mutation rate p <-> P
		) {

    // Extra calculations
    const double pref_level_het = pow(pref_level, dominance);

    const double no_t_mut_rate = 1.0 - t_mut_rate;
    const double no_p_mut_rate = 1.0 - p_mut_rate;

    const double prob_xxx_mut = no_p_mut_rate * no_t_mut_rate * no_t_mut_rate;
    const double prob_pxx_mut =    p_mut_rate * no_t_mut_rate * no_t_mut_rate;
    const double prob_xtx_mut = no_p_mut_rate *    t_mut_rate * no_t_mut_rate;
    const double prob_ptx_mut =    p_mut_rate *    t_mut_rate * no_t_mut_rate;
    //const double prob_xxt_mut = no_p_mut_rate * no_t_mut_rate *    t_mut_rate;  // Duplicate of xtx
    //const double prob_pxt_mut =    p_mut_rate * no_t_mut_rate *    t_mut_rate;  // Duplicate of ptx
    const double prob_xtt_mut = no_p_mut_rate *    t_mut_rate *    t_mut_rate;
    const double prob_ptt_mut =    p_mut_rate *    t_mut_rate *    t_mut_rate;

    // Determne the fitness of each genotype
    const double f_TT_fit = 1.0 + female_gain;
    const double f_tT_fit = 1.0 + dominance * female_gain;
    const double f_tt_fit = 1.0;
    const double m_TT_fit = 1.0 - male_harm;
    const double m_tT_fit = 1.0 - dominance * male_harm;
    const double m_tt_fit = 1.0;

    pop_values next_pop;
    Normalize();
    for (size_t UD = 0; UD < updates; UD++) {
      // Determine the weight of each female genotype.
      double f_pTT_weight = female_pTT * f_TT_fit;
      double f_ptT_weight = female_ptT * f_tT_fit;
      double f_ptt_weight = female_ptt * f_tt_fit;  // Technically, just * 1.0
      double f_PTT_weight = female_PTT * f_TT_fit;
      double f_PtT_weight = female_PtT * f_tT_fit;
      double f_Ptt_weight = female_Ptt * f_tt_fit;  // Technically, just * 1.0

      // Determine the weight of particular gene combinations being inherited from the female.
      double f_Pt_weight = f_PtT_weight / 2.0 + f_Ptt_weight;            // a
      double f_P_weight  = f_PTT_weight + f_PtT_weight + f_Ptt_weight;   // b
      double f_PT_weight = f_PTT_weight + f_PtT_weight / 2.0;            // c
      double f_pt_weight = f_ptT_weight / 2.0 + f_ptt_weight;            // d
      double f_p_weight  = f_pTT_weight + f_ptT_weight + f_ptt_weight;   // e
      double f_pT_weight = f_pTT_weight + f_ptT_weight / 2.0;            // g

      // Determine the weight of each male genotype
      double m_TT_weight = male_TT * m_TT_fit;
      double m_tT_weight = male_tT * m_tT_fit;
      double m_tt_weight = male_tt * m_tt_fit;

      // Determine how preferences help...
      double ave_pref_level = pref_level * male_TT + pref_level_het * male_tT + male_tt;
      double TT_pref_ratio = pref_level / ave_pref_level;
      double tT_pref_ratio = pref_level_het / ave_pref_level;
      double tt_pref_ratio = 1.0 / ave_pref_level;

      // Determine weight of each female in the next population.
      double next_f_pTT = f_pT_weight * m_TT_weight + 
                          f_pT_weight * m_tT_weight / 2.0;
      double next_f_ptT = f_pt_weight * m_TT_weight + 
                          f_p_weight  * m_tT_weight / 2.0 +
                          f_pT_weight * m_tt_weight;
      double next_f_ptt = f_pt_weight * m_tT_weight / 2.0 + 
                          f_pt_weight * m_tt_weight;
      double next_f_PTT = f_PT_weight * m_TT_weight * TT_pref_ratio + 
                          f_PT_weight * m_tT_weight * tT_pref_ratio / 2.0;
      double next_f_PtT = f_Pt_weight * m_TT_weight * TT_pref_ratio +
                          f_P_weight  * m_tT_weight * tT_pref_ratio / 2.0 +
                          f_PT_weight * m_tt_weight * tt_pref_ratio;
      double next_f_Ptt = f_Pt_weight * m_tT_weight * tT_pref_ratio / 2.0 +
                          f_Pt_weight * m_tt_weight * tt_pref_ratio;

      // Now account for flow due to mutations!
      next_pop.female_pTT = next_f_pTT * prob_xxx_mut + next_f_PTT * prob_pxx_mut +
                            next_f_ptT * prob_xtx_mut + next_f_PtT * prob_ptx_mut + 
                            next_f_ptt * prob_xtt_mut + next_f_Ptt * prob_ptt_mut;

      next_pop.female_ptT = (next_f_pTT * prob_xtx_mut + next_f_PTT * prob_ptx_mut) * 2.0 +  // Either t can toggle
                            next_f_ptT * prob_xxx_mut + next_f_PtT * prob_pxx_mut + 
                            (next_f_ptt * prob_xtx_mut + next_f_Ptt * prob_ptx_mut) * 2.0 +
                            next_f_ptT * prob_xtt_mut + next_f_PtT * prob_ptt_mut;           // BOTH t's can toggle

      next_pop.female_ptt = next_f_pTT * prob_xtt_mut + next_f_PTT * prob_ptt_mut +
                            next_f_ptT * prob_xtx_mut + next_f_PtT * prob_ptx_mut + 
                            next_f_ptt * prob_xxx_mut + next_f_Ptt * prob_pxx_mut;

      next_pop.female_PTT = next_f_pTT * prob_pxx_mut + next_f_PTT * prob_xxx_mut +
                            next_f_ptT * prob_ptx_mut + next_f_PtT * prob_xtx_mut + 
                            next_f_ptt * prob_ptt_mut + next_f_Ptt * prob_xtt_mut;

      next_pop.female_PtT = (next_f_pTT * prob_ptx_mut + next_f_PTT * prob_xtx_mut) * 2.0 +
                            next_f_ptT * prob_pxx_mut + next_f_PtT * prob_xxx_mut + 
                            (next_f_ptt * prob_ptx_mut + next_f_Ptt * prob_xtx_mut) * 2.0 +
                            next_f_ptT * prob_ptt_mut + next_f_PtT * prob_xtt_mut;

      next_pop.female_Ptt = next_f_pTT * prob_ptt_mut + next_f_PTT * prob_xtt_mut +
                            next_f_ptT * prob_ptx_mut + next_f_PtT * prob_xtx_mut + 
                            next_f_ptt * prob_pxx_mut + next_f_Ptt * prob_xxx_mut;

      // Clean up results.
      next_pop.Normalize();
      *this = next_pop;
    }
  }
};


int main() {
  pop_values pop;

  for (size_t i = 0; i < 5000000; i++) {    
    if (i % 1000 == 0) pop.MinPrint(i);
    pop.Run();
  }
  std::cout << "Final:\n";
  pop.Print();
}
