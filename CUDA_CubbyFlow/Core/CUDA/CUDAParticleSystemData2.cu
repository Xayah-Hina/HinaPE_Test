// This code is based on Jet framework.
// Copyright (c) 2018 Doyub Kim
// CubbyFlow is voxel-based fluid simulation engine for computer games.
// Copyright (c) 2020 CubbyFlow Team
// Core Part: Chris Ohk, Junwoo Hwang, Jihong Sin, Seungwoo Yoo
// AI Part: Dongheon Cho, Minseo Kim
// We are making my contributions/submissions to this project solely in our
// personal capacity and are not conveying any rights to any intellectual
// property of any third parties.

#include "CUDAParticleSystemData2Func.hpp"

#include <Core/CUDA/CUDAParticleSystemData2.hpp>
#include <Core/CUDA/ThrustUtils.hpp>

#include <thrust/copy.h>
#include <thrust/scan.h>
#include <thrust/transform.h>

using namespace CubbyFlow;

constexpr uint32_t DEFAULT_HASH_GRID_RESOLUTION = 64;

CUDAParticleSystemData2::CUDAParticleSystemData2() : CUDAParticleSystemData2(0)
{
    // Do nothing
}

CUDAParticleSystemData2::CUDAParticleSystemData2(size_t numberOfParticles)
{
    m_positionIdx = AddVectorData();
    m_velocityIdx = AddVectorData();

    Resize(numberOfParticles);
}

CUDAParticleSystemData2::CUDAParticleSystemData2(
    const CUDAParticleSystemData2& other)
{
    Set(other);
}

CUDAParticleSystemData2::CUDAParticleSystemData2(
    CUDAParticleSystemData2&& other) noexcept
    : m_numberOfParticles(std::exchange(other.m_numberOfParticles, 0)),
      m_positionIdx(std::exchange(other.m_positionIdx, 0)),
      m_velocityIdx(std::exchange(other.m_velocityIdx, 0)),
      m_intDataList(std::move(other.m_intDataList)),
      m_floatDataList(std::move(other.m_floatDataList)),
      m_vectorDataList(std::move(other.m_vectorDataList)),
      m_neighborSearcher(std::move(other.m_neighborSearcher)),
      m_neighborStarts(std::move(other.m_neighborStarts)),
      m_neighborEnds(std::move(other.m_neighborEnds)),
      m_neighborLists(std::move(other.m_neighborLists))
{
    // Do nothing
}

CUDAParticleSystemData2& CUDAParticleSystemData2::operator=(
    const CUDAParticleSystemData2& other)
{
    Set(other);
    return *this;
}

CUDAParticleSystemData2& CUDAParticleSystemData2::operator=(
    CUDAParticleSystemData2&& other) noexcept
{
    m_numberOfParticles = std::exchange(other.m_numberOfParticles, 0);
    m_positionIdx = std::exchange(other.m_positionIdx, 0);
    m_velocityIdx = std::exchange(other.m_velocityIdx, 0);

    m_intDataList = std::move(other.m_intDataList);
    m_floatDataList = std::move(other.m_floatDataList);
    m_vectorDataList = std::move(other.m_vectorDataList);

    m_neighborSearcher = std::move(other.m_neighborSearcher);
    m_neighborStarts = std::move(other.m_neighborStarts);
    m_neighborEnds = std::move(other.m_neighborEnds);
    m_neighborLists = std::move(other.m_neighborLists);
    return *this;
}

void CUDAParticleSystemData2::Resize(size_t newNumberOfParticles)
{
    m_numberOfParticles = newNumberOfParticles;

    for (auto& attr : m_intDataList)
    {
        attr.Resize(newNumberOfParticles, 0);
    }

    for (auto& attr : m_floatDataList)
    {
        attr.Resize(newNumberOfParticles, 0.0f);
    }

    for (auto& attr : m_vectorDataList)
    {
        attr.Resize(newNumberOfParticles, make_float2(0, 0));
    }
}

size_t CUDAParticleSystemData2::NumberOfParticles() const
{
    return m_numberOfParticles;
}

size_t CUDAParticleSystemData2::AddIntData(int initialVal)
{
    const size_t attrIdx = m_intDataList.size();
    m_intDataList.emplace_back(NumberOfParticles(), initialVal);
    return attrIdx;
}

size_t CUDAParticleSystemData2::AddFloatData(float initialVal)
{
    const size_t attrIdx = m_floatDataList.size();
    m_floatDataList.emplace_back(NumberOfParticles(), initialVal);
    return attrIdx;
}

size_t CUDAParticleSystemData2::AddVectorData(const Vector2F& initialVal)
{
    const size_t attrIdx = m_vectorDataList.size();
    m_vectorDataList.emplace_back(NumberOfParticles(), ToFloat2(initialVal));
    return attrIdx;
}

size_t CUDAParticleSystemData2::NumberOfIntData() const
{
    return m_intDataList.size();
}

size_t CUDAParticleSystemData2::NumberOfFloatData() const
{
    return m_floatDataList.size();
}

size_t CUDAParticleSystemData2::NumberOfVectorData() const
{
    return m_vectorDataList.size();
}

ConstCUDAArrayView1<float2> CUDAParticleSystemData2::Positions() const
{
    return VectorDataAt(m_positionIdx);
}

CUDAArrayView1<float2> CUDAParticleSystemData2::Positions()
{
    return VectorDataAt(m_positionIdx);
}

ConstCUDAArrayView1<float2> CUDAParticleSystemData2::Velocities() const
{
    return VectorDataAt(m_velocityIdx);
}

CUDAArrayView1<float2> CUDAParticleSystemData2::Velocities()
{
    return VectorDataAt(m_velocityIdx);
}

ConstCUDAArrayView1<int> CUDAParticleSystemData2::IntDataAt(size_t idx) const
{
    return m_intDataList[idx].View();
}

CUDAArrayView1<int> CUDAParticleSystemData2::IntDataAt(size_t idx)
{
    return m_intDataList[idx].View();
}

ConstCUDAArrayView1<float> CUDAParticleSystemData2::FloatDataAt(
    size_t idx) const
{
    return m_floatDataList[idx].View();
}

CUDAArrayView1<float> CUDAParticleSystemData2::FloatDataAt(size_t idx)
{
    return m_floatDataList[idx].View();
}

ConstCUDAArrayView1<float2> CUDAParticleSystemData2::VectorDataAt(
    size_t idx) const
{
    return m_vectorDataList[idx].View();
}

CUDAArrayView1<float2> CUDAParticleSystemData2::VectorDataAt(size_t idx)
{
    return m_vectorDataList[idx].View();
}

void CUDAParticleSystemData2::AddParticle(const Vector2F& newPosition,
                                          const Vector2F& newVelocity)
{
    std::vector<float2> hostPos;
    std::vector<float2> hostVel;
    hostPos.push_back(ToFloat2(newPosition));
    hostVel.push_back(ToFloat2(newVelocity));

    CUDAArray1<float2> devicePos{ hostPos };
    CUDAArray1<float2> deviceVel{ hostVel };

    AddParticles(devicePos, deviceVel);
}

void CUDAParticleSystemData2::AddParticles(
    ConstArrayView1<Vector2F> newPositions,
    ConstArrayView1<Vector2F> newVelocities)
{
    std::vector<float2> hostPos(newPositions.Length());
    std::vector<float2> hostVel(newVelocities.Length());

    for (size_t i = 0; i < newPositions.Length(); ++i)
    {
        hostPos[i] = ToFloat2(newPositions[i]);
    }
    for (size_t i = 0; i < newVelocities.Length(); ++i)
    {
        hostVel[i] = ToFloat2(newVelocities[i]);
    }

    CUDAArray1<float2> devicePos{ hostPos };
    CUDAArray1<float2> deviceVel{ hostVel };

    AddParticles(devicePos, deviceVel);
}

void CUDAParticleSystemData2::AddParticles(
    ConstCUDAArrayView1<float2> newPositions,
    ConstCUDAArrayView1<float2> newVelocities)
{
    assert(newVelocities.Length() == newPositions.Length());

    size_t oldNumberOfParticles = NumberOfParticles();
    Resize(oldNumberOfParticles + newPositions.Length());

    auto pos = Positions();
    CUDACopy(newPositions.data(), newPositions.Length(),
             pos.data() + oldNumberOfParticles);

    if (newVelocities.Length() > 0)
    {
        auto vel = Velocities();
        CUDACopy(newVelocities.data(), newVelocities.Length(),
                 vel.data() + oldNumberOfParticles);
    }
}

ConstCUDAArrayView1<uint32_t> CUDAParticleSystemData2::NeighborStarts() const
{
    return m_neighborStarts.View();
}

ConstCUDAArrayView1<uint32_t> CUDAParticleSystemData2::NeighborEnds() const
{
    return m_neighborEnds.View();
}

const CUDAPointHashGridSearcher2* CUDAParticleSystemData2::NeighborSearcher()
    const
{
    return m_neighborSearcher.get();
}

ConstCUDAArrayView1<uint32_t> CUDAParticleSystemData2::NeighborLists() const
{
    return m_neighborLists.View();
}

void CUDAParticleSystemData2::BuildNeighborSearcher(float maxSearchRadius)
{
    if (m_neighborSearcher == nullptr)
    {
        m_neighborSearcher = std::make_shared<CUDAPointHashGridSearcher2>(
            DEFAULT_HASH_GRID_RESOLUTION, DEFAULT_HASH_GRID_RESOLUTION,
            2.0f * maxSearchRadius);
    }

    m_neighborSearcher->Build(Positions());
}

void CUDAParticleSystemData2::BuildNeighborLists(float maxSearchRadius)
{
    m_neighborStarts.Resize(m_numberOfParticles);
    m_neighborEnds.Resize(m_numberOfParticles);

    auto neighborStarts = m_neighborStarts.View();

    // Count nearby points
    thrust::for_each(
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(0) + NumberOfParticles(),
        ForEachNeighborFunc<NoOpFunc, CountNearbyPointsFunc>(
            *m_neighborSearcher, maxSearchRadius, Positions().data(),
            NoOpFunc(), CountNearbyPointsFunc(m_neighborStarts.data())));

    // Make start/end point of neighbor list, and allocate neighbor list.
    thrust::inclusive_scan(thrustCBegin(m_neighborStarts),
                           thrustCEnd(m_neighborStarts),
                           thrustBegin(m_neighborEnds));
    thrust::transform(thrustCBegin(m_neighborEnds), thrustCEnd(m_neighborEnds),
                      thrustCBegin(m_neighborStarts),
                      thrustBegin(m_neighborStarts),
                      thrust::minus<unsigned int>());
    size_t rbeginIdx =
        m_neighborEnds.Length() > 0 ? m_neighborEnds.Length() - 1 : 0;
    uint32_t m = m_neighborEnds[rbeginIdx];
    m_neighborLists.Resize(m, 0);

    // Build neighbor lists
    thrust::for_each(
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(0) + NumberOfParticles(),
        ForEachNeighborFunc<BuildNeighborListsFunc, NoOpFunc>(
            *m_neighborSearcher, maxSearchRadius, Positions().data(),
            BuildNeighborListsFunc(m_neighborStarts.data(),
                                   m_neighborEnds.data(),
                                   m_neighborLists.data()),
            NoOpFunc()));
}

void CUDAParticleSystemData2::Set(const CUDAParticleSystemData2& other)
{
    m_numberOfParticles = other.m_numberOfParticles;
    m_positionIdx = other.m_positionIdx;
    m_velocityIdx = other.m_velocityIdx;

    m_intDataList = other.m_intDataList;
    m_floatDataList = other.m_floatDataList;
    m_vectorDataList = other.m_vectorDataList;

    if (other.m_neighborSearcher != nullptr)
    {
        m_neighborSearcher = std::make_shared<CUDAPointHashGridSearcher2>(
            *other.m_neighborSearcher);
    }

    m_neighborStarts = other.m_neighborStarts;
    m_neighborEnds = other.m_neighborEnds;
    m_neighborLists = other.m_neighborLists;
}
