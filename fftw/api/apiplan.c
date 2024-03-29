/*
 * Copyright (c) 2003, 2007 Matteo Frigo
 * Copyright (c) 2003, 2007 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "api.h"

static plan *mkplan0(planner *plnr, unsigned flags, 
		     const problem *prb, int hash_info, 
		     wisdom_state_t wisdom_state)
WITH_ALIGNED_STACK({
     /* map API flags into FFTW flags */
     X(mapflags)(plnr, flags);

     plnr->flags.hash_info = hash_info;
     plnr->wisdom_state = wisdom_state;

     /* create plan */
     return plnr->adt->mkplan(plnr, prb);
})

static void aligned_awake(plan *ego, enum wakefulness wakefulness)
WITH_ALIGNED_STACK({
     X(plan_awake)(ego, wakefulness);
})

static unsigned force_estimator(unsigned flags)
{
     flags &= ~(FFTW_MEASURE | FFTW_PATIENT | FFTW_EXHAUSTIVE);
     return (flags | FFTW_ESTIMATE);
}

static plan *mkplan(planner *plnr, unsigned flags, 
		    const problem *prb, int hash_info)
{
     plan *pln;

     pln = mkplan0(plnr, flags, prb, hash_info, WISDOM_NORMAL);

     if (plnr->wisdom_state == WISDOM_NORMAL && !pln) {
	  /* maybe the planner failed because of inconsistent wisdom;
	     plan again ignoring infeasible wisdom */
	  pln = mkplan0(plnr, force_estimator(flags), prb, 
			hash_info, WISDOM_IGNORE_INFEASIBLE);
     }

     if (plnr->wisdom_state == WISDOM_IS_BOGUS) {
	  /* if the planner detected a wisdom inconsistency,
	     forget all wisdom and plan again */
	  plnr->adt->forget(plnr, FORGET_EVERYTHING);

	  A(!pln);
	  pln = mkplan0(plnr, flags, prb, hash_info, WISDOM_NORMAL);

	  if (plnr->wisdom_state == WISDOM_IS_BOGUS) {
	       /* if it still fails, plan without wisdom */
	       plnr->adt->forget(plnr, FORGET_EVERYTHING);

	       A(!pln);
	       pln = mkplan0(plnr, force_estimator(flags), 
			     prb, hash_info, WISDOM_IGNORE_ALL);
	  }
     }

     return pln;
}

apiplan *X(mkapiplan)(int sign, unsigned flags, problem *prb)
{
     apiplan *p = 0;
     plan *pln;
     unsigned flags_used_for_planning;
     planner *plnr = X(the_planner)();
     unsigned int pats[] = {FFTW_ESTIMATE, FFTW_MEASURE,
			    FFTW_PATIENT, FFTW_EXHAUSTIVE};
     int pat, pat_max;

     if (flags & FFTW_WISDOM_ONLY) {
	  /* special mode that returns a plan only if wisdom is present,
	     and returns 0 otherwise.  Phil Dumont requested this
	     feature, do not remove. */
	  flags_used_for_planning = flags;
	  pln = mkplan0(plnr, flags, prb, 0, WISDOM_ONLY);
     } else {
	  pat_max = flags & FFTW_ESTIMATE ? 0 :
	       (flags & FFTW_EXHAUSTIVE ? 3 :
		(flags & FFTW_PATIENT ? 2 : 1));
	  pat = plnr->timelimit >= 0 ? 0 : pat_max;

	  flags &= ~(FFTW_ESTIMATE | FFTW_MEASURE | 
		     FFTW_PATIENT | FFTW_EXHAUSTIVE);

	  plnr->start_time = X(get_crude_time)();
	  
	  /* plan at incrementally increasing patience until we run
	     out of time */
	  for (pln = 0, flags_used_for_planning = 0; pat <= pat_max; ++pat) {
	       plan *pln1;
	       unsigned tmpflags = flags | pats[pat];
	       pln1 = mkplan(plnr, tmpflags, prb, 0);

	       if (!pln1) {
		    /* don't bother continuing if planner failed or timed out */
		    A(!pln || plnr->timed_out);
		    break;
	       }

	       X(plan_destroy_internal)(pln);
	       pln = pln1;
	       flags_used_for_planning = tmpflags;
	  }
     }

     if (pln) {
	  /* build apiplan */
	  p = (apiplan *) MALLOC(sizeof(apiplan), PLANS);
	  p->prb = prb;
	  p->sign = sign; /* cache for execute_dft */
	  
	  /* re-create plan from wisdom, adding blessing */
	  p->pln = mkplan(plnr, flags_used_for_planning, prb, BLESSING);

	  if (sizeof(trigreal) > sizeof(R)) {
	       /* this is probably faster, and we have enough trigreal
		  bits to maintain accuracy */
	       aligned_awake(p->pln, AWAKE_SQRTN_TABLE);
	  } else {
	       /* more accurate */
	       aligned_awake(p->pln, AWAKE_SINCOS);
	  }
	  
	  /* we don't use pln for p->pln, above, since by re-creating the
	     plan we might use more patient wisdom from a timed-out mkplan */
	  X(plan_destroy_internal)(pln);
     } else
	  X(problem_destroy)(prb);
     
     /* discard all information not necessary to reconstruct the plan */
     plnr->adt->forget(plnr, FORGET_ACCURSED);
     
     return p;
}

void X(destroy_plan)(X(plan) p)
{
     if (p) {
          X(plan_awake)(p->pln, SLEEPY);
          X(plan_destroy_internal)(p->pln);
          X(problem_destroy)(p->prb);
          X(ifree)(p);
     }
}
