#include <gmpxx.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <chrono>
#include <sstream>
using Q=mpq_class;
struct Poly {
 std::vector<Q> a;
 Poly(){} Poly(long n){if(n)a.push_back(Q(n));} Poly(const Q& n){if(n!=0)a.push_back(n);}
 void trim(){while(!a.empty() && a.back()==0)a.pop_back();}
 bool zero()const{return a.empty();} bool one()const{return a.size()==1 && a[0]==1;}
 int deg()const{return (int)a.size()-1;}
 Poly monic()const{if(zero())return *this;Poly p=*this;Q q=a.back();for(auto&x:p.a)x/=q;return p;}
 std::string str()const{std::string s="[";for(size_t i=0;i<a.size();i++){if(i)s+=",";s+="\""+a[i].get_str()+"\"";}return s+"]";}
};
bool operator==(const Poly&a,const Poly&b){return a.a==b.a;}
Poly operator+(const Poly&a,const Poly&b){Poly c;c.a.resize(std::max(a.a.size(),b.a.size()));for(size_t i=0;i<a.a.size();i++)c.a[i]+=a.a[i];for(size_t i=0;i<b.a.size();i++)c.a[i]+=b.a[i];c.trim();return c;}
Poly operator-(const Poly&a,const Poly&b){Poly c;c.a.resize(std::max(a.a.size(),b.a.size()));for(size_t i=0;i<a.a.size();i++)c.a[i]+=a.a[i];for(size_t i=0;i<b.a.size();i++)c.a[i]-=b.a[i];c.trim();return c;}
Poly operator*(const Poly&a,const Poly&b){if(a.zero()||b.zero())return Poly();Poly c;c.a.resize(a.a.size()+b.a.size()-1);for(size_t i=0;i<a.a.size();i++)for(size_t j=0;j<b.a.size();j++)c.a[i+j]+=a.a[i]*b.a[j];c.trim();return c;}
Poly scale(Poly a,const Q&s){for(auto&x:a.a)x*=s;a.trim();return a;}
std::pair<Poly,Poly> divrem(Poly a,const Poly&b){if(b.zero())throw std::runtime_error("divide zero");Poly q;if(a.deg()<b.deg())return {q,a};q.a.resize(a.deg()-b.deg()+1);while(!a.zero()&&a.deg()>=b.deg()){int k=a.deg()-b.deg();Q c=a.a.back()/b.a.back();q.a[k]+=c;for(int j=0;j<=b.deg();j++)a.a[k+j]-=c*b.a[j];a.trim();}q.trim();return {q,a};}
Poly exact_div(const Poly&a,const Poly&b){auto qr=divrem(a,b);if(!qr.second.zero())throw std::runtime_error("nonexact division");return qr.first;}
Poly gcd(Poly a,Poly b){if(a.deg()==0||b.deg()==0)return Poly(1);while(!b.zero()){Poly r=divrem(a,b).second;a=b;b=r;}return a.monic();}
Poly MOD;
Poly reduce_mod(const Poly&a){return MOD.zero()?a:divrem(a,MOD).second;}
Poly inverse_mod(Poly a){a=reduce_mod(a);Poly b=MOD,s0(1),s1(0);while(!b.zero()){auto qr=divrem(a,b);a=b;b=qr.second;Poly sn=s0-qr.first*s1;s0=s1;s1=sn;}if(a.deg()!=0)throw std::runtime_error("nonunit in residue field");return reduce_mod(scale(s0,1/a.a[0]));}
struct RF {
 Poly n,d;
 RF():n(),d(1){} RF(long a):n(a),d(1){} RF(const Poly&a):n(a),d(1){if(!MOD.zero())n=reduce_mod(n);}
 RF(Poly a,Poly b):n(a),d(b){norm();}
 bool zero()const{return n.zero();}
 void norm(){if(!MOD.zero()){n=reduce_mod(n);d=reduce_mod(d);if(d.zero())throw std::runtime_error("residue divzero");if(!d.one())n=reduce_mod(n*inverse_mod(d));d=Poly(1);return;}if(n.zero()){d=Poly(1);return;}if(d.zero())throw std::runtime_error("RF divzero");if(d.deg()==0){n=scale(n,1/d.a[0]);d=Poly(1);return;}Poly g=gcd(n,d);if(!g.one()){n=exact_div(n,g);d=exact_div(d,g);}Q c=d.a.back();if(c!=1){n=scale(n,1/c);d=scale(d,1/c);}}
 RF inverse()const{return RF(d,n);}
};
RF operator-(const RF&a){RF b=a;b.n=scale(b.n,-1);return b;}
RF operator+(const RF&a,const RF&b){if(a.zero())return b;if(b.zero())return a;if(a.d.one()&&b.d.one())return RF(a.n+b.n);if(a.d==b.d)return RF(a.n+b.n,a.d);return RF(a.n*b.d+b.n*a.d,a.d*b.d);}
RF operator-(const RF&a,const RF&b){if(b.zero())return a;return a+(-b);}
RF operator*(const RF&a,const RF&b){if(a.zero()||b.zero())return RF();if(a.d.one()&&b.d.one())return RF(a.n*b.n);return RF(a.n*b.n,a.d*b.d);}
struct MaskCmp {bool operator()(int a,int b)const{int ca=__builtin_popcount((unsigned)a),cb=__builtin_popcount((unsigned)b);return ca<cb||(ca==cb&&a<b);}};
using Vec=std::map<int,RF,MaskCmp>;
Vec mul(const Vec&a,const Vec&b){Vec c;for(const auto&[i,x]:a)for(const auto&[j,y]:b)if(!(i&j)){int k=i|j;RF z=x*y;auto it=c.find(k);if(it==c.end()){if(!z.zero())c[k]=z;}else{z=it->second+z;if(z.zero())c.erase(it);else it->second=z;}}return c;}
struct Space {
 std::map<int,Vec> rows;std::set<std::string> bad;
 void track(const RF&a){if(!MOD.zero())return;if(a.n.deg()>0)bad.insert(a.n.monic().str());if(a.d.deg()>0)bad.insert(a.d.monic().str());}
 bool insert(Vec &v){while(!v.empty()){int p=v.begin()->first;RF c=v.begin()->second;auto it=rows.find(p);if(it==rows.end()){track(c);RF inv=c.inverse();for(auto&[k,x]:v)x=x*inv;rows[p]=v;return true;}for(const auto&[k,x]:it->second){RF z=c*x;auto jt=v.find(k);if(jt==v.end()){if(!z.zero())v[k]=-z;}else{z=jt->second-z;if(z.zero())v.erase(jt);else jt->second=z;}}}return false;}
};
bool ZERO_QUAD=false;
std::vector<Vec> generators(const std::vector<std::pair<int,int>>&edges){std::vector<Vec> gg;Poly zz;zz.a={0,1};RF b(zz);for(int v=1;v<=8;v++){std::vector<std::pair<int,int>>inc;for(size_t i=0;i<edges.size();i++){auto e=edges[i];if(e.first==v||e.second==v)inc.push_back({1<<i,v==std::min(e.first,e.second)?1:-1});}Vec g;for(int sub=1;sub<(1<<(int)inc.size());sub++){int k=__builtin_popcount((unsigned)sub),m=0,s=1;for(size_t j=0;j<inc.size();j++)if(sub>>j&1){m|=inc[j].first;s*=inc[j].second;}if(k==1)g[m]=RF(s);if(k==2&&!ZERO_QUAD)g[m]=RF(2*s);if(k==3)g[m]=RF(6*s)*b;}gg.push_back(g);}return gg;}
int main(int argc,char**argv){int which=argc>1?std::stoi(argv[1]):1;
if(argc>2){std::stringstream ss(argv[2]);std::string term;while(std::getline(ss,term,',')){Q a(term);a.canonicalize();MOD.a.push_back(a);}MOD.trim();}
if(argc>3)ZERO_QUAD=std::string(argv[3])=="zero";
std::vector<std::pair<int,int>>edges={{1,2},{2,3},{3,4},{4,1},{1,5},{5,6},{2,7},{which==1?3:4,8}};auto gg=generators(edges);Space sp;Vec unit;unit[0]=RF(1);sp.insert(unit);std::vector<Vec>front={unit};std::vector<int>hh={1};for(int d=1;d<=8;d++){auto start=std::chrono::steady_clock::now();std::vector<Vec>next;for(const auto&v:front)for(const auto&g:gg){Vec w=mul(v,g);if(sp.insert(w))next.push_back(w);}if(next.empty())break;hh.push_back(next.size());std::cout<<"degree "<<d<<" h "<<next.size()<<" total "<<sp.rows.size()<<" seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;front=next;if(sp.rows.size()==240)break;}
std::cout<<"H [";for(size_t i=0;i<hh.size();i++){if(i)std::cout<<",";std::cout<<hh[i];}std::cout<<"]"<<std::endl;std::cout<<"PIVOTS [";bool first=true;for(const auto&s:sp.bad){if(!first)std::cout<<",";first=false;std::cout<<s;}std::cout<<"]"<<std::endl;return 0;}
