#include <bits/stdc++.h>

void debug_out() { std::cerr << std::endl; }
template <typename Head, typename... Tail>
void debug_out(Head H, Tail... T) {
    std::cerr << " " << H;
    debug_out(T...);
}
#ifdef LOCAL
#define dbg(...) std::cerr << "[" << #__VA_ARGS__ << "]:", debug_out(__VA_ARGS__)
#else
#define dbg(...)
#endif
#define all(x) (x).begin(), (x).end()
#define sz(x) ((int)(x).size())
#define pii std::pair<int, int>
using namespace std;
using i64 = int64_t;


int INF = INT_MAX - 67, MIN_SEED = 10;
vector<int> ks = {3, 6, 9};
string fa_path = "GCF_000005845.2_ASM584v2_genomic.fna", fq_path = "ERR022075_1.fastq";

int n;
string t;
vector<int> sa;


struct hit {
    int ed, l, r, rc, loci;
};



vector<int> SUFFIX_ARRAY(string& s) 
{

    int n = sz(s), A = 128;
    vector<int> p(n), c(n), cnt(max(A, n), 0), pn(n), cn(n);

    for (int i = 0; i < n; i++) {
        cnt[s[i]]++;
    }

    for (int i = 1; i < A; i++) {
        cnt[i] += cnt[i - 1];
    }

    for (int i = n - 1; i >= 0; i--) {
        p[--cnt[s[i]]] = i;
    }

    int cls = 1;
    c[p[0]] = 0;
    for (int i = 1; i < n; i++) {
        if (s[p[i]] != s[p[i - 1]]) {
            cls++;
        }

        c[p[i]] = cls - 1;

    }

    for (int h = 1; h < n && cls < n; h <<= 1) {

        for (int i = 0; i < n; i++) {
            pn[i] = p[i] - h;
            if (pn[i] < 0) {
                pn[i] += n;
            }
        }

        fill(cnt.begin(), cnt.begin() + cls, 0);

        for (int i = 0; i < n; i++) {
            cnt[c[pn[i]]]++;
        }

        for (int i = 1; i < cls; i++) {
            cnt[i] += cnt[i - 1];
        }

        for (int i = n - 1; i >= 0; i--) {
            p[--cnt[c[pn[i]]]] = pn[i];
        }

        cn[p[0]] = 0;
        cls = 1;

        for (int i = 1; i < n; i++) {
            int a = p[i], b = p[i - 1];
            if (c[a] != c[b] || c[(a + h) % n] != c[(b + h) % n]) {
                ++cls;
            }

            cn[a] = cls - 1;
        }

        swap(c, cn);
    }

    return p;
}


string REV(string& s) {
    string r(s.rbegin(), s.rend());

    for (auto& c: r) {
        if (c == 'A') {
            c = 'T';
        } 
        else if (c == 'T') {
            c = 'A';
        } 
        else if (c == 'C') {
            c = 'G';
        } 
        else if (c == 'G') {
            c = 'C';
        }
    }
    return r;
}


void find(string& r, int k, vector<int>& c) 
{
    int len = sz(r) / (k + 1);


    for (int q = 0; q <= k; q++) {
        int o = q * len, lo = 0, hi = n + 1;

        while (lo < hi) {

            int mid = (lo + hi) / 2;

            if (t.compare(sa[mid], len, r, o, len) < 0) {
                lo = mid + 1;
            } 
            
            else {
                hi = mid;
            }
        }

        while (lo <= n && t.compare(sa[lo], len, r, o, len) == 0) {
            c.push_back(sa[lo] - o);
            lo++;
        }

    }


}


hit check(string& r, int s, int k) {
    int m = sz(r), w = 2 * k + 3;
    int dp[2][64], st[2][64];

    fill(dp[0], dp[0]+w, INF);
    fill(dp[1], dp[1]+w, INF);

    for (int b = 1; b <= 2 * k + 1; b++) {
        dp[0][b] = 0;
        st[0][b] = s - k + b-1;
    }


    for (int i = 0; i < m; i++) 
    {
        int *pd = dp[i&1], *ps = st[i&1], *nd = dp[(i+1) & 1], *ns = st[(i+1) & 1];
        int mn = INF;

        for (int b = 1; b <= 2 * k + 1; b++) 
        {
            int j = s + i - k + b-1;
            nd[b] = pd[b] + (j >= 0 && j < n && t[j] == r[i]? 0: 1);
            ns[b] = ps[b];
            if (pd[b + 1] + 1 < nd[b]) {
                nd[b] = pd[b + 1] + 1;
                ns[b] = ps[b + 1];
            }

            if (nd[b - 1] + 1 < nd[b]) {
                nd[b] = nd[b - 1] + 1;
                ns[b] = ns[b - 1];
            }

            mn = min(mn, nd[b]);
        }


        if (mn > k) {
            return {k+1, 0, 0, 0, 0};
        }

    }

    int *fd = dp[m&1], *fs = st[m & 1], bb = 1;
    for (int b = 1; b <= 2 * k + 1; b++) {
        if (fd[b] < fd[bb]) {
            bb = b;
        }
    }

    return {fd[bb], fs[bb], s+m-k + bb-1, 0, 0};

}



hit base(string& r, string& rc, int k, i64 id) {
    vector<hit> best;
    vector<int> cand;


    for (int z = 0; z < 2; z++) 
    {
        string& q = (z? rc: r);
        cand.clear();

        find(q, k, cand);
        sort(all(cand));

        cand.erase(unique(all(cand)), cand.end());

        int lst = -INF;
        for (int s: cand) 
        {
            hit cur = check(q, s, k);
            if (cur.ed > k) {
                continue;
            }

            cur.rc = z;

            if (best.empty() || cur.ed < best[0].ed) {
                best.clear();
                best.push_back(cur);
                lst = cur.r;
            } 

            else if (cur.ed == best[0].ed && cur.l >= lst) {
                best.push_back(cur);
                lst = cur.r;
            }

        }
    
    }

    if (best.empty()) {
        return {INF, 0, 0, 0, 0};
    }

    hit res = best[id%sz(best)];

    res.loci = sz(best);
    return res;

}



hit align(string& r, i64 id) {
    string rc = REV(r);
    for (int k: ks) {
        if (sz(r) / (k + 1) < MIN_SEED) {
            break;
        }
        hit res = base(r, rc, k, id);
        if (res.ed <= k) {
            return res;
        }
    }

    return {INF, 0, 0, 0, 0};
}


void solve() {

    ifstream fa(fa_path), fq(fq_path);

    string line;
    while (getline(fa, line)) 
    {
        if (line.empty() || line[0] == '>') {
            continue;
        }

        for (char ch: line) {
            if (isalpha(ch)) {
                t += toupper(ch);
            }
        }
    }

    n = sz(t);
    t += '$';
    sa = SUFFIX_ARRAY(t);

    double t_idx = 1.0 * clock() / CLOCKS_PER_SEC;
    cerr << "genome " << n << " bp, SA built in " << t_idx << " s\n";

    int kmx = ks.back();
    i64 tot = 0, mapped = 0, uniq = 0, multi = 0, rev = 0, len_sum = 0, ed_sum = 0, q_map = 0, q_unm = 0, b_map = 0, b_unm = 0;
    vector<i64> cnt_ed(kmx+1, 0);
    vector<int> d_all(n+1, 0), d_uniq(n+1, 0);
    string name, r, plus, qual;



    while (getline(fq, name) && getline(fq, r) && getline(fq, plus) && getline(fq, qual)) 
    {
        hit h = align(r, tot);
        tot++;
        int m = sz(r);
        i64 q_sum = 0;
        for (char ch: qual) {
            q_sum += ch - 33;
        }

        if (h.ed <= kmx) {
            mapped++;
            cnt_ed[h.ed]++;
            ed_sum += h.ed;
            len_sum += m;
            rev += h.rc;
            q_map += q_sum;
            b_map += m;
            int l = max(h.l, 0), rr = min(h.r, n);
            d_all[l]++;
            d_all[rr]--;

            if (h.loci == 1) {
                uniq++;
                d_uniq[l]++;
                d_uniq[rr]--;
            } 

            else {
                multi++;
            }
        }

        else {
            q_unm += q_sum;
            b_unm += m;
        }

        if (tot % 1000000 == 0) {
            cerr << tot << " reads, mapped " << 100.0 * mapped / tot << "%\n";
        }

    }





    double t_all = 1.0 * clock() / CLOCKS_PER_SEC;



    auto pct = [](i64 a, i64 b) { return b? 100.0 * a / b: 0.0; };
    cout << fixed << setprecision(2);
    cout << "genome length: " << n << "\n";
    cout << "threshold: edit distance <= " << kmx << " over the whole read\n";
    cout << "total reads: " << tot << "\n";
    cout << "mapped reads: " << mapped << " (" << pct(mapped, tot) << "%)\n";
    cout << "unmapped reads: " << tot - mapped << " (" << pct(tot - mapped, tot) << "%)\n";
    cout << "unique: " << uniq << " (" << pct(uniq, tot) << "% of all, " << pct(uniq, mapped) << "% of mapped)\n";
    cout << "multi-mapped: " << multi << " (" << pct(multi, tot) << "% of all, " << pct(multi, mapped) << "% of mapped)\n";
    cout << "forward / reverse: " << mapped - rev << " / " << rev << "\n";

    cout << "\nedit distance distribution (mapped reads):\n";
    i64 acc = 0;
    for (int e = 0; e <= kmx; e++) {
        acc += cnt_ed[e];
        cout << "ed = " << e << ": " << cnt_ed[e] << " (" << pct(cnt_ed[e], tot) << "% of all, cumulative " << pct(acc, tot) << "%)\n";
    }
    cout << "mean edit distance: " << (mapped? 1.0 * ed_sum / mapped: 0.0) << "\n";
    cout << "mean identity: " << (len_sum? 100.0 * (1.0 - 1.0 * ed_sum / len_sum): 0.0) << "%\n";
    cout << "mean base quality: mapped " << (b_map? 1.0 * q_map / b_map: 0.0) << ", unmapped " << (b_unm? 1.0 * q_unm / b_unm: 0.0) << "\n";




    for (int z = 0; z < 2; z++) {
        vector<int>& d = (z? d_uniq: d_all);
        i64 sum = 0, c1 = 0, c10 = 0;
        int cur = 0, mx = 0;
        vector<int> dep(n);
        for (int i = 0; i < n; i++) {
            cur += d[i];
            dep[i] = cur;
            sum += cur;
            c1 += cur >= 1;
            c10 += cur >= 10;
            mx = max(mx, cur);
        }
        nth_element(dep.begin(), dep.begin() + n / 2, dep.end());
        cout << "\ncoverage (" << (z? "unique reads only": "all mapped reads") << "):\n";
        cout << "covered bases (>= 1x): " << c1 << " / " << n << " (" << pct(c1, n) << "%)\n";
        cout << "covered bases (>= 10x): " << c10 << " (" << pct(c10, n) << "%)\n";
        cout << "mean depth: " << 1.0 * sum / n << "x, median " << dep[n / 2] << "x, max " << mx << "x\n";
    }

    cout << "\nindex time: " << t_idx << " s\n";
    cout << "total time: " << t_all << " s\n";


}



int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);


    int t = 1;
    // cin >> t;
    while (t--) {
        solve();
    }

    return 0;
}