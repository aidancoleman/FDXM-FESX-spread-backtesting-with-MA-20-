// csv_two_cols.c
// compiled as :  gcc -O2 -std=c11 csv_two_cols.c -o cols -lm
// Ran as :      ./cols "./demo "<path to csv>" --skip 21" to skip 21 first rows of data as we use MA(20)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <errno.h>
#include <math.h>


/* percentile helper calculate running 5th and 95th percentile values for 100 sample window (sliding window starting from end of array and decrementing backwards back to start of array) */

static void reverse_copy(const double *src, double *dst, size_t n){
    for (size_t i = 0; i < n; ++i) dst[i] = src[n - 1 - i];
}

static int cmp_double_inc(const void *a, const void *b){
    const double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y);
}

//forward trade simulator -> no overlaps in trades
//X = our series of prices we are looking at
//dir = direction of trade, +1 -> long, -1 -> short
//stop_offset = stop size in units of X (get out at stop if breached)
//take_profit = target
// min_hold (time of minimum hold) and stop/Target are in the same units as X !
static long simulate_trade_forward(const double *X, size_t n, size_t k, int dir, size_t min_hold, double stop_offset, double take_profit, double *out_exit_price, double *out_profit)
{
    //if entry is out of bounds then abort
    if (k >= n || !isfinite(X[k])) return -1;
    const double entry = X[k];
    bool reached_target = false; // becomes true once cur crosses the target in our favour
    
    double stop_long  = entry - stop_offset; // for dir=+1, -> long means we stop if price gets lower
    double stop_short = entry + stop_offset; // for dir=-1, -> short means we stop if price gets higher than this
    
    size_t t = k; //index of last filled bar
    size_t held = 0; //how many bars since entry (for minimum hold time)
    //iterate forward one bar at a time (don't peek, simulates actual trading with no future knowledge)
    while (t + 1 < n) {
        size_t tnext = t + 1;
        while (tnext < n && !isfinite(X[tnext])) tnext++;
        if (tnext >= n) break;

        double cur = X[tnext];
        if (dir > 0) {
            if (cur >= take_profit) reached_target = true;
        } else {
            if (cur <= take_profit) reached_target = true;
        }
        
        // hard stop
        if (dir > 0) {
            if (cur <= stop_long) { *out_exit_price = cur; *out_profit = (cur - entry); return (long)tnext; }
        } else {
            if (cur >= stop_short){ *out_exit_price = cur; *out_profit = (entry - cur); return (long)tnext; }
        }

        // move to breakeven after +1R (stop size), so once unrealised profit exceeds our stop in the other direction we raise our stop to the entry price
        // new: lock in +0.5R once price has moved +1R
        static const double R_BE = 1.0;   // trigger at +1R
        static const double R_LOCK = 0.5; // move stop to +0.5R

        if (dir > 0 && cur >= entry + R_BE*stop_offset) {
            // long: stop below current price; raise it above entry by 0.5R
            stop_long = entry + R_LOCK * stop_offset;
        } else if (dir < 0 && cur <= entry - R_BE*stop_offset) {
            // short: stop above current price; lower it below entry by 0.5R
            stop_short = entry - R_LOCK * stop_offset;
        }

        /*fixed take-profit (requires min_hold satisfied)
        if (held + 1 >= min_hold) {
            if (dir > 0) {
                if (cur >= take_profit) { *out_exit_price = cur; *out_profit = (cur - entry); return (long)tnext; }
            } else {
                if (cur <= take_profit) { *out_exit_price = cur; *out_profit = (entry - cur); return (long)tnext; }
            }
        }*/

        // fallback exit after min_hold on first adverse tick, so if it ticks back against us after holding for the minimum hold time (if winning) we will exit
        if (held + 1 >= min_hold && reached_target) {
            double prev = X[t];
            if (isfinite(prev)) {
                if (dir > 0 && cur < prev)  { *out_exit_price = cur; *out_profit = (cur - entry); return (long)tnext; }
                if (dir < 0 && cur > prev)  { *out_exit_price = cur; *out_profit = (entry - cur); return (long)tnext; }
            }
        }

        held++;
        t = tnext;
    }

    if (isfinite(X[t])) {
        *out_exit_price = X[t];
        *out_profit = (dir > 0) ? (X[t] - entry) : (entry - X[t]);
        return (long)t;
    }
    return -1;
}

// x_sorted must be sorted ascending, m >= 1, p in [0,1]
static double percentile_inc_from_sorted(const double *x_sorted, size_t m, double p){
    if (m == 0) return NAN;
    if (m == 1) return x_sorted[0];
    if (p <= 0) return x_sorted[0];
    if (p >= 1) return x_sorted[m-1];
    double r = p * (double)(m - 1);         // 0-based fractional index
    size_t lo = (size_t)floor(r);
    size_t hi = (size_t)ceil(r);
    double a = x_sorted[lo], b = x_sorted[hi];
    double frac = r - (double)lo;
    return a + frac * (b - a);
}
//converting excel columns to correct csv indexing
static int excel_col_to_index(const char *s){
    long idx=0; for(const char *p=s; *p; ++p){ if(!isalpha((unsigned char)*p)) return -1; idx=idx*26+(toupper((unsigned char)*p)-'A'+1); }
    return (int)(idx-1);
}
static char **read_csv_record(FILE *fp, int *out_count){
    *out_count=0; size_t cap=0; char *line=NULL; ssize_t nread=getline(&line,&cap,fp);
    if(nread<0){ free(line); return NULL; }
    size_t fcap=16, flen=0; char **fields=malloc(fcap*sizeof *fields); if(!fields){ free(line); return NULL; }
    const char *p=line;
    while(*p){
        if(flen==fcap){ fcap*=2; char **nf=realloc(fields, fcap*sizeof *nf); if(!nf){ goto fail; } fields=nf; }
        size_t bcap=64, blen=0; char *buf=malloc(bcap); if(!buf) goto fail;
        bool quoted=false; if(*p=='"'){ quoted=true; p++; }
        for(;;){
            char c=*p++;
            if(quoted){
                if(c=='"'){ if(*p=='"'){ if(blen+1>=bcap){ bcap*=2; buf=realloc(buf,bcap); if(!buf) goto fail; } buf[blen++]='"'; p++; }
                            else quoted=false; }
                else if(c=='\0') break;
                else { if(blen+1>=bcap){ bcap*=2; buf=realloc(buf,bcap); if(!buf) goto fail; } buf[blen++]=c; }
            }else{
                if(c==','||c=='\n'||c=='\r'||c=='\0'){ if(c=='\r' && *p=='\n') p++; break; }
                else { if(blen+1>=bcap){ bcap*=2; buf=realloc(buf,bcap); if(!buf) goto fail; } buf[blen++]=c; }
            }
        }
        buf[blen]='\0'; fields[flen++]=buf;
        if(*p==',') p++; while(*p=='\r'||*p=='\n') p++; if(*p=='\0') break;
    }
    free(line); *out_count=(int)flen; return fields;
fail:
    if(fields){ for(size_t i=0;i<flen;++i) free(fields[i]); free(fields); }
    free(line); return NULL;
}
//failure check for over/under flow
static double parse_or_nan(const char *s){
    //null and blank guard
    if(!s) return NAN; while(*s && isspace((unsigned char)*s)) s++; if(*s=='\0') return NAN;
    errno=0; char *endp=NULL; double x=strtod(s,&endp); if(errno!=0 || endp==s) return NAN;
    //endp points to first unparsed character
    while(*endp && isspace((unsigned char)*endp)) endp++; return *endp? NAN : x;
}


//statement to execute on trigger, print the index and current value of dist. from MA(20) and spread values
static void print_trigger(size_t i, double Bi, double Di, double p5, double p95, void *user){
    (void)user;
    printf("[TRIGGER] i=%zu  B[i]=%g  D[i]=%g  p5(Bwin)=%g  p95(Bwin)=%g\n", i, Bi, Di, p5, p95);
}

/*reusable loader*/
int csv_load_two_cols(const char *path, const char *col1_spec, const char *col2_spec, long skip_lines, double **out1, double **out2, size_t *out_n)
{
    int c1 = isalpha((unsigned char)col1_spec[0]) ? excel_col_to_index(col1_spec) : atoi(col1_spec);
    int c2 = isalpha((unsigned char)col2_spec[0]) ? excel_col_to_index(col2_spec) : atoi(col2_spec);
    if(c1<0 || c2<0) return -2;

    FILE *fp=fopen(path,"r"); if(!fp) return -1;

    // skip top lines
    for(long i=0;i<skip_lines;++i){
        int nf=0; char **f=read_csv_record(fp,&nf); if(!f) break;
        for(int j=0;j<nf;++j) free(f[j]); free(f);
    }

    size_t cap=1024, n=0;
    double *a1=malloc(cap*sizeof *a1), *a2=malloc(cap*sizeof *a2);
    if(!a1||!a2){ fclose(fp); free(a1); free(a2); return -3; }

    int nf=0; char **fields=NULL;
    while((fields=read_csv_record(fp,&nf))!=NULL){
        if (n == cap) {
            cap *= 2;
            double *t1 = realloc(a1, cap * sizeof *t1);
            if (!t1) { perror("realloc a1"); break; }
            double *t2 = realloc(a2, cap * sizeof *t2);
            if (!t2) { perror("realloc a2"); a1 = t1; break; }
            a1 = t1; a2 = t2;
        }
        a1[n] = (c1<nf)? parse_or_nan(fields[c1]) : NAN;
        a2[n] = (c2<nf)? parse_or_nan(fields[c2]) : NAN;
        n++;
        for(int i=0;i<nf;++i) free(fields[i]); free(fields);
    }
    fclose(fp);

    *out1=a1; *out2=a2; *out_n=n;
    return 0;
}

// Trade record struct
typedef struct {
    size_t entry_idx, exit_idx;
    int dir;                // +1 long, -1 short
    double entry_price, exit_price;
    double profit;          // dir * (exit - entry)
    double p5_at_entry, p95_at_entry;
} Trade;

// Forward trade simulator on D, starting at k.
// Uses min_hold (e.g., 20) and stop_offset (in D units).
// Returns exit index (>=k), fills *out_profit and *out_exit_price, or -1 if no viable data ahead.
static long simulate_trade_forward_D(const double *D, size_t n, size_t k, int dir, size_t min_hold, double stop_offset, double take_profit, double *out_exit_price, double *out_profit)
{
    if (k >= n || !isfinite(D[k])) return -1;
    const double entry = D[k];
    bool reached_target = false; // becomes true once cur crosses the target in our favor
    // Stop levels in D units
    double stop_long  = entry - stop_offset; // for dir=+1
    double stop_short = entry + stop_offset; // for dir=-1

    size_t t = k;
    size_t held = 0;

    // advance one step at a time (no peeking)
    while (t + 1 < n) {
        // advance to next finite
        size_t tnext = t + 1;
        while (tnext < n && !isfinite(D[tnext])) tnext++;
        if (tnext >= n) break;

        // check stop at the new price
        double cur = D[tnext];
        /*
        if (dir > 0) {
            if (cur >= take_profit) reached_target = true;
        } else {
            if (cur <= take_profit) reached_target = true;
        }*/
        
        if (dir > 0) { // long
            if (cur <= stop_long) {
                *out_exit_price = cur;
                *out_profit = (cur - entry);
                return (long)tnext;
            }
        } else {// short
            if (cur >= stop_short) {
                *out_exit_price = cur;
                *out_profit = (entry - cur);
                return (long)tnext;
            }
        }
        
        // new: lock in +0.5R once price has moved +1R
        static const double R_BE = 1.0;   // trigger at +1R
        static const double R_LOCK = 0.5; // move stop to +0.5R

        if (dir > 0 && cur >= entry + R_BE*stop_offset) {
            // long: stop below current price; raise it above entry by 0.5R
            stop_long = entry + R_LOCK * stop_offset;
        } else if (dir < 0 && cur <= entry - R_BE*stop_offset) {
            // short: stop above current price; lower it below entry by 0.5R
            stop_short = entry - R_LOCK * stop_offset;
        }
        
        if (dir > 0) {
            if (cur >= take_profit && held + 1 >= min_hold) {
                *out_exit_price = cur; *out_profit = (cur - entry);
                return (long)tnext;
            }
        } else {
            if (cur <= take_profit && held + 1 >= min_hold) {
                *out_exit_price = cur; *out_profit = (entry - cur);
                return (long)tnext;
            }
        }

        // after min_hold, exit on first adverse tick (guard against reversal)
        if (held + 1 >= min_hold && reached_target) {
            // compare with previous finite (D[t])
            double prev = D[t];
            if (isfinite(prev)) {
                if (dir > 0) { // long: exit on first down-tick
                    if (cur < prev) {
                        *out_exit_price = cur;
                        *out_profit = (cur - entry);
                        return (long)tnext;
                    }
                } else {        // short: exit on first up-tick
                    if (cur > prev) {
                        *out_exit_price = cur;
                        *out_profit = (entry - cur);
                        return (long)tnext;
                    }
                }
            }
        }

        // keep holding
        held++;
        t = tnext;
    }

    // If we run out of data, exit at last seen finite
    if (isfinite(D[t])) {
        *out_exit_price = D[t];
        *out_profit = (dir > 0) ? (D[t] - entry) : (entry - D[t]);
        return (long)t;
    }
    return -1;
}

// Forward scan with non-overlapping trades.
// At time t (>= window-1), compute p5/p95 of B[t-window+1...t]. If B[t] breaches:
//   B[t] <= p5  -> LONG D at t
//   B[t] >= p95 -> SHORT D at t
// Then simulate forward on D until exit; resume scanning at the bar AFTER the exit.
// Forward scan, NON-overlapping.
// Percentiles computed on P (D), trades executed on T (B).
size_t run_strategy_forward_nonoverlap(
    const double *P,          // percentile series (D)
    const double *T,          // trade series (B)
    size_t n,
    size_t window, double p_low, double p_high,
    size_t min_hold,
    Trade **out_trades, size_t *out_m
){
    *out_trades = NULL; *out_m = 0;
    if (!P || !T || n == 0 || window == 0 || n < window) return 0;

    double *buf = malloc(window * sizeof *buf);
    if (!buf) { perror("malloc"); return 0; }

    size_t cap = 16, m = 0;
    Trade *trades = malloc(cap * sizeof *trades);
    if (!trades) { free(buf); perror("malloc"); return 0; }

    for (size_t t = window - 1; t < n; ++t) {
        // window on P = D
        size_t start = t - (window - 1);
        size_t w = 0;
        for (size_t k = start; k <= t; ++k)
            if (isfinite(P[k])) buf[w++] = P[k];

        if (w == 0 || !isfinite(P[t]) || !isfinite(T[t])) continue;

        qsort(buf, w, sizeof *buf, cmp_double_inc);
        double p5  = percentile_inc_from_sorted(buf, w, p_low);
        double p95 = percentile_inc_from_sorted(buf, w, p_high);

        int dir = 0;
        // strict breach; change to <=/>= for inclusive of percentile values
        if (isfinite(p5)  && P[t] <  p5) dir = +1;      // LONG on T(B)
        else if (isfinite(p95) && P[t] >  p95) dir = -1; // SHORT on T(B)

        if (dir != 0) {
            // Stops/Target are in T(B) units
            double S = 10.0; // fixed stop in spread units
            double R = 2.0; // TPrice multiple
            double stop_offset = S;
            double take_profit = (dir > 0) ? (T[t] + R*S) : (T[t] - R*S);

            double exit_px = NAN, profit = NAN;
            long jexit = simulate_trade_forward(T, n, t, dir, min_hold, stop_offset, take_profit, &exit_px, &profit);
            if (jexit < 0) break;

            if (m == cap) {
                cap *= 2;
                Trade *tmp = realloc(trades, cap * sizeof *tmp);
                if (!tmp) { perror("realloc"); break; }
                trades = tmp;
            }
            trades[m++] = (Trade){
                .entry_idx = t,
                .exit_idx  = (size_t)jexit,
                .dir = dir,
                .entry_price = T[t],     // spread pts at entry
                .exit_price  = exit_px,  // spread pts at exit
                .profit      = (dir > 0 ? (exit_px - T[t]) : (T[t] - exit_px)),
                .p5_at_entry = p5,
                .p95_at_entry = p95
            };

            // non-overlap: jump to exit bar
            t = (size_t)jexit;
        }
    }

    free(buf);
    *out_trades = trades; *out_m = m;
    return m;
}
// Backward scan: 100-window percentiles of B; on breach at i, trade D starting at i.
size_t run_strategy_backward(
    const double *B, const double *D, size_t n,
    size_t window, // 100 for us
    double p_low, double p_high, // 0.05 and 0.95
    size_t min_hold, // e.g., 20
    Trade **out_trades, size_t *out_m // returns heap array of trades
){
    *out_trades = NULL; *out_m = 0;
    if (!B || !D || n < window) return 0;

    // temp buffer for spread pts from MA(20) window
    double *buf = malloc(window * sizeof *buf);
    if (!buf) { perror("malloc"); return 0; }

    size_t cap = 16, m = 0;
    Trade *trades = malloc(cap * sizeof *trades);
    if (!trades) { free(buf); perror("malloc"); return 0; }

    // scan from end down to window-1
    for (size_t i = n - 1; i + 1 >= window; ) {
        size_t start = i - (window - 1);

        // collect finite B values in the window
        size_t w = 0;
        for (size_t k = start; k <= i; ++k) if (isfinite(B[k])) buf[w++] = B[k];
        if (w == 0 || !isfinite(B[i]) || !isfinite(D[i])) {
            if (i == 0) break; else { --i; continue; }
        }

        qsort(buf, w, sizeof *buf, cmp_double_inc);
        double p5 = percentile_inc_from_sorted(buf, w, p_low);
        double p95 = percentile_inc_from_sorted(buf, w, p_high);

        int dir = 0;
        if (isfinite(p5) && B[i] <= p5) dir = +1; // go long D
        else if (isfinite(p95) && B[i] >= p95) dir = -1; // go short D

        if (dir != 0) {
            //choose stop size (in spread pts units)
        
            double S = 25.0, R = 2.0;
            double stop_offset = S;
            double take_profit = (dir > 0) ? (D[i] + R*S) : (D[i] - R*S);

            // simulate forward trade on spread from i
            double exit_px = NAN, profit = NAN;
            long jexit = simulate_trade_forward_D(D, n, i, dir, min_hold, stop_offset, take_profit, &exit_px, &profit);
            if (jexit >= 0) {
                if (m == cap) { cap *= 2; Trade *t2 = realloc(trades, cap * sizeof *t2); if (!t2) { perror("realloc"); break; } trades = t2; }
                trades[m++] = (Trade){
                    .entry_idx = i,
                    .exit_idx = (size_t)jexit,
                    .dir = dir,
                    .entry_price = D[i],
                    .exit_price  = exit_px,
                    .profit = (dir > 0 ? (exit_px - D[i]) : (D[i] - exit_px)),
                    .p5_at_entry = p5, .p95_at_entry = p95
                };
                // resume scanning at the bar before the entry to avoid overlapping detections
                if (i == 0) break; else --i;
            } else {
                // no exit found; still step back
                if (i == 0) break; else --i;
            }
        } else {
            if (i == 0) break; else --i;
        }
    }

    free(buf);
    *out_trades = trades; *out_m = m;
    return m;
}


/* load spread & distance from MA(20) after skipping 20 lines (blanks), then carry out algo with sliding window */
static size_t count_finite(const double *x, size_t n){ size_t k=0; for(size_t i=0;i<n;++i) if(isfinite(x[i])) ++k; return k; }
static double mean_finite(const double *x, size_t n){ double s=0; size_t k=0; for(size_t i=0;i<n;++i) if(isfinite(x[i])){ s+=x[i]; ++k; } return k? s/k : NAN; }

int main(int argc, char **argv){
    if (argc < 2){
        fprintf(stderr,"Usage: %s <csv_path> [--skip N]\n", argv[0]);
        return 1;
    }
    const char *path = argv[1];
    long skip = 0;
    if (argc >= 4 && strcmp(argv[2], "--skip") == 0) {
        skip = strtol(argv[3], NULL, 10);
        if (skip < 0) skip = 0;
    } else if (argc == 3 && strncmp(argv[2], "--skip=", 7) == 0) {
        skip = strtol(argv[2] + 7, NULL, 10);
        if (skip < 0) skip = 0;
    }

    double *colB=NULL, *colD=NULL; size_t n=0;
    int rc = csv_load_two_cols(path, "B", "D", skip, &colB, &colD, &n);
    if (rc){ fprintf(stderr,"load failed (rc=%d)\n", rc); return 1; }

    // Strategy parameters
    size_t W = 100;          // backward window
    double p_low  = 0.05;
    double p_high = 0.95;
    size_t min_hold = 20;

    // reverse arrays so index increases with time
    double *Brev = malloc(n * sizeof *Brev);
    double *Drev = malloc(n * sizeof *Drev);
    if (!Brev || !Drev) { perror("malloc rev"); free(Brev); free(Drev); free(colB); free(colD); return 1; }
    reverse_copy(colB, Brev, n);
    reverse_copy(colD, Drev, n);

    // Run on chronological (reversed) arrays
    Trade *trades_rev = NULL; size_t m = 0;
    run_strategy_forward_nonoverlap(Drev, Brev, n, W, p_low, p_high, min_hold, &trades_rev, &m);

    // Report, mapping indices back to original array indexing!!
    double total = 0.0;
    for (size_t t = 0; t < m; ++t) {
        const Trade *tr = &trades_rev[t];
        size_t entry_orig = n - 1 - tr->entry_idx;  // map back
        size_t exit_orig  = n - 1 - tr->exit_idx;   // map back

        printf("Trade %zu: %s  entry(orig i=%zu) D=%.6g  exit(orig j=%zu) D=%.6g  profit=%.6g  (p5=%.6g, p95=%.6g)\n",
               t+1, tr->dir>0?"LONG":"SHORT",
               entry_orig, colD[entry_orig],
               exit_orig,  colD[exit_orig],
               tr->profit, tr->p5_at_entry, tr->p95_at_entry);

        total += tr->profit;
    }
    printf("Total trades: %zu, total P&L: %.6g\n", m, total);

    // cleanup
    free(trades_rev);
    free(Brev); free(Drev);
    free(colB); free(colD);
    return 0;
}
