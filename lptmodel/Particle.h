#ifndef PARTICLE_H_
#define PARTICLE_H_
#include <Eigen/Dense>
#include "io.h"
#include "macros.h"
#include "Fluid.h"
#include "typedefs.h"
#include <memory>
#include "DynamicFactory.h"
#include "ParticleForces.h"

/* Particle base class */
class Particle {
public:
	virtual ~Particle() { }
	virtual Particle * clone() const = 0;
	Particle() = default;
	Particle(const Particle &);
	Particle & operator=(const Particle &);

	// Getters and setters
	GETSET(Vector, position)
	GETSET(Vector, velocity)
	GETSET(Matrix, shear)
	GETSET(scalar, pas)
	GETSET(scalar, dose)
	GETSET(scalar, age)
	GETSET(bool, isAlive)
	GETSET(int, id)
	GETSET(int, collisionCount)
	GETSET(scalar, injectionTime)
	GETSET(Vector, fluidVelocity)
	GETSET(Vector, fluidVorticity)

	GETSET(scalar, density)
	GETSET(scalar, radius)


	GETSET(std::vector<std::unique_ptr<ParticleForce>>, particleForces)

	virtual void updateMomentum(scalar dt, const Vector & fluidVelocity, const Matrix & shear, const Fluid & fluid, const Vector & fluidVorticity = {0., 0., 0.}) = 0;
	virtual int typeId() const = 0;
	virtual std::string typeName() const = 0;
	virtual void fromJSON(const json & jsonObject) { }
	virtual void writeBinary(std::ostream & os) const;
	virtual void readBinary(std::istream & is);

private:
	Vector position_{0., 0., 0.};
	Vector velocity_{0., 0., 0.};
	int id_{-1};
	Matrix shear_;
	scalar pas_{0.};
	scalar dose_{0.};
	scalar age_{0.};
	bool isAlive_{true};
	int collisionCount_{0};
	scalar injectionTime_{-1};
	Vector fluidVelocity_{0., 0., 0.};
	Vector fluidVorticity_{0., 0., 0.};
	std::vector<std::unique_ptr<ParticleForce>> particleForces_;

protected:
	scalar density_{0};
	scalar radius_{0};
};

/* Particle creation */
using ParticleFactory = DynamicFactory<Particle>;

ParticleFactory & particleFactory();

template<class DerivedType>
inline int registerParticleTypeToFactory(const std::string & typeName)
{
	return particleFactory().registerCreator(new NamedSimpleObjectCreator<Particle, DerivedType>(typeName));
}

/* Tracer particle */
class TracerParticle : public Particle {
public:
	TracerParticle * clone() const override;
	void updateMomentum(scalar dt, const Vector & fluidVelocity, const Matrix & shear, const Fluid & fluid, const Vector & fluidVorticity = {0., 0., 0.}) override;
	int typeId() const override { return typeId_; }
	std::string typeName() const override { return typeName_; }

private:
	static int typeId_;
	static std::string typeName_;
};

/* Material particle */
class MaterialParticle : public Particle {
public:
	// MaterialParticle() = default;
	// MaterialParticle(const MaterialParticle &);
	// MaterialParticle & operator=(const MaterialParticle &);


	
	

	MaterialParticle * clone() const override;
	void updateMomentum(scalar dt, const Vector & fluidVelocity, const Matrix & shear, const Fluid & fluid, const Vector & fluidVorticity = {0., 0., 0.}) override;
	int typeId() const override;
	std::string typeName() const override { return typeName_;} 
	void fromJSON(const json & jsonObject) override;
	void writeBinary(std::ostream & out) const override;
	void readBinary(std::istream & in) override;
	

private:
	static int typeId_;
	static std::string typeName_;
};

/* Particle moving with constant velocity (for testing) */
class NoMomentumUpdateParticle : public Particle {
public:
	NoMomentumUpdateParticle * clone() const { return new NoMomentumUpdateParticle(*this); }
	void updateMomentum(scalar dt, const Vector & fluidVelocity, const Matrix & shear, const Fluid & fluid, const Vector & fluidVorticity = {0., 0., 0.}) { }
	int typeId() const override { return typeId_; }
	std::string typeName() const override { return typeName_;} 
	

private:
	static int typeId_;
	static std::string typeName_;
};


#endif /* PARTICLE_H_ */
