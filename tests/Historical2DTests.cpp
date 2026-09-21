#include "solvers/mac2d/ReflectionMacFluidSolver2D.hpp"
#include "fixtures/Ba20987Smoke2D.hpp"
#include <fstream>
#include <iostream>
#include <string>

void snapshot(const std::vector<float>& smoke, int n, const std::string& path) {
    const int margin=std::min(32,(n-1)/2), size=n-2*margin;
    std::ofstream out(path,std::ios::binary);
    out<<"P6\n"<<size<<' '<<size<<"\n255\n";
    for(int y=n-margin-1;y>=margin;--y) for(int x=margin;x<n-margin;++x) {
        const float a=std::clamp(smoke[x+n*y]*2.5f,0.f,1.f);
        for (auto channel : {std::pair{.03f,.78f},std::pair{.04f,.80f},std::pair{.07f,.87f}})
            out.put(static_cast<char>(std::clamp(int(255*(channel.first+(channel.second-channel.first)*a)),0,255)));
    }
    if(!out) throw std::runtime_error("Snapshot write failed");
}

int main(int argc,char** argv) {
    try {
        // A quarter-cell translation has a known linear-interpolation result.
        // This ensures SL actually bypasses the MacCormack correction.
        ReflectionMacFluidSolver2D sl(16), mc(16);
        sl.macCormackSmoke = false;
        for (auto* solver : {&sl, &mc}) {
            std::fill(solver->vx.begin(), solver->vx.end(), 1.f);
            solver->smoke[solver->idC(7, 8)] = 1.f;
            solver->advectScalar(solver->smoke, .25f * solver->h());
        }
        for (int y=0;y<16;++y) for (int x=0;x<16;++x) {
            const float expected = y==8 ? (x==7 ? .75f : x==8 ? .25f : 0.f) : 0.f;
            if (std::abs(sl.smoke[sl.idC(x,y)]-expected)>1e-6f)
                throw std::runtime_error("Reflection SL translation failed");
        }
        if (sl.smoke == mc.smoke) throw std::runtime_error("Advection selection has no effect");
        sl.resetState();
        for (int f=0;f<120;++f) {
            sl.advance(1.f/60.f);
            for (float d : sl.smoke) if (!std::isfinite(d) || d<0.f || d>1.f)
                throw std::runtime_error("Reflection SL smoke bounds");
            for (int i=0;i<sl.n;++i)
                if (sl.vx[sl.idX(0,i)]!=0.f || sl.vx[sl.idX(sl.n,i)]!=0.f ||
                    sl.vy[sl.idY(i,0)]!=0.f || sl.vy[sl.idY(i,sl.n)]!=0.f)
                    throw std::runtime_error("Reflection SL wall leakage");
        }
        const auto slSmoke = sl.smoke, slVx = sl.vx, slVy = sl.vy;
        sl.resetState();
        if (sl.macCormackSmoke) throw std::runtime_error("Reset changed advection selection");
        for (int f=0;f<120;++f) sl.advance(1.f/60.f);
        if (sl.smoke!=slSmoke || sl.vx!=slVx || sl.vy!=slVy)
            throw std::runtime_error("Reflection SL replay differs");
        const int n=argc>1?std::stoi(argv[1]):24;
        const int frames=argc>2?std::stoi(argv[2]):24;
        ReflectionMacFluidSolver2D current(n);
        // No-penetration with free tangential motion, not no-slip zero velocity.
        std::fill(current.vx.begin(),current.vx.end(),.7f);
        std::fill(current.vy.begin(),current.vy.end(),-.3f);
        current.applyBoundary();
        for(int i=0;i<n;++i) {
            if(current.vx[current.idX(0,i)]!=0.f || current.vx[current.idX(n,i)]!=0.f ||
               current.vy[current.idY(i,0)]!=0.f || current.vy[current.idY(i,n)]!=0.f)
                throw std::runtime_error("Normal wall velocity is not zero");
        }
        if(current.vx[current.idX(n/2,0)]!=.7f || current.vx[current.idX(n/2,n-1)]!=.7f ||
           current.vy[current.idY(0,n/2)]!=-.3f || current.vy[current.idY(n-1,n/2)]!=-.3f)
            throw std::runtime_error("Wall condition suppressed tangential slip");
        current.resetState();
        ba20987::SmokeSim2D original(n);
        original.stirStrength=0.f;
        original.openLeft=original.openRight=original.openTop=original.openBottom=false;
        double maximumDifference=0.;
        for(int frame=1;frame<=frames;++frame) {
            original.step(1.f/60.f,.8f,.995f,1.f,60);
            current.advance(1.f/60.f);
            for(int i=0;i<n;++i)
                if(current.vx[current.idX(0,i)]!=0.f || current.vx[current.idX(n,i)]!=0.f ||
                   current.vy[current.idY(i,0)]!=0.f || current.vy[current.idY(i,n)]!=0.f)
                    throw std::runtime_error("Reflection step leaked through a wall");
            for(auto pair : {std::pair{&original.smoke,&current.smoke},std::pair{&original.vx,&current.vx},
                            std::pair{&original.vy,&current.vy},std::pair{&original.pressure,&current.pressure}}) {
                for(std::size_t i=0;i<pair.first->size();++i) {
                    const float a=(*pair.first)[i], b=(*pair.second)[i];
                    if(!std::isfinite(a) || !std::isfinite(b)) throw std::runtime_error("Non-finite historical field");
                    maximumDifference=std::max(maximumDifference,double(std::abs(a-b)));
                }
            }
            if(maximumDifference!=0.) throw std::runtime_error("Historical numerical parity failed: "+std::to_string(maximumDifference));
            if(frame%60==0 || frame==frames) {
                std::cout<<"n="<<n<<" frame="<<frame<<" max_field_difference="<<maximumDifference
                         <<" divergence="<<current.diagnostics().maxDivergence<<std::endl;
                if(argc>3) {
                    snapshot(original.smoke,n,std::string(argv[3])+"-old-"+std::to_string(frame)+".ppm");
                    snapshot(current.smoke,n,std::string(argv[3])+"-new-"+std::to_string(frame)+".ppm");
                }
            }
        }
        const auto expected=current.smoke;
        current.resetState();
        // Short replay also verifies reset and deterministic advance.
        if(frames<=24) {
            for(int f=0;f<frames;++f) current.advance(1.f/60.f);
            if(current.smoke!=expected) throw std::runtime_error("Reset replay differs");
        }
        std::cout<<"PASS reflection parity; closed free-slip walls and stirring disabled in both\n";
    } catch(const std::exception& e) { std::cerr<<"FAIL "<<e.what()<<'\n'; return 1; }
}
