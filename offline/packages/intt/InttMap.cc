#include "InttMap.h"

#include <boost/format.hpp>

InttMap::Online_s& InttMap::Online_s::operator++()
{
	if(!InttMap::IsValidOrWildcard(*this))
	{
		*this = InttMap::OnlineEnd;
		return *this;
	}

	if(chn != InttMap::Wildcard)
	{
		if(++chn < 128)
		{
			return *this;
		}
		chn = 0;
	}
	if(chp != InttMap::Wildcard)
	{
		if(++chp < 26)
		{
			return *this;
		}
		chp = 0;
	}
	if(arm != InttMap::Wildcard)
	{
		if(++arm < 2)
		{
			return *this;
		}
		arm = 0;
	}
	if(ldr != InttMap::Wildcard)
	{
		if(++ldr < (lyr < 2 ? 12 : 16))
		{
			return *this;
		}
		ldr = 0;
	}
	if(lyr != InttMap::Wildcard)
	{
		if(++lyr < 4)
		{
			return *this;
		}
		lyr = 0;
	}

	*this = InttMap::OnlineEnd;
	return *this;
}

InttMap::Online_s& InttMap::Online_s::operator--()
{
	if(!InttMap::IsValidOrWildcard(*this))
	{
		*this = InttMap::OnlineREnd;
		return *this;
	}

	if(chn != InttMap::Wildcard)
	{
		if(--chn >= 0)
		{
			return *this;
		}
		chn = 127;
	}
	if(chp != InttMap::Wildcard)
	{
		if(--chp >= 0)
		{
			return *this;
		}
		chp = 25;
	}
	if(arm != InttMap::Wildcard)
	{
		if(--arm >= 0)
		{
			return *this;
		}
		arm = 1;
	}
	if(ldr != InttMap::Wildcard)
	{
		if(--ldr >= 0)
		{
			return *this;
		}
		ldr = lyr < 3 ? 11 : 15;
	}
	if(lyr != InttMap::Wildcard)
	{
		if(--lyr >= 0)
		{
			return *this;
		}
		lyr = 3;
	}

	*this = InttMap::OnlineREnd;
	return *this;
}

bool InttMap::IsValid(
	InttMap::Online_s const& onl)
{
	if(onl.lyr < 0 || 3 < onl.lyr)
	{
		return false;
	}
	if(onl.ldr < 0 || (onl.lyr < 2 ? 11 : 15) < onl.ldr)
	{
		return false;
	}
	if(onl.arm < 0 || 1 < onl.arm)
	{
		return false;
	}
	if(onl.chp < 0 || 25 < onl.chp)
	{
		return false;
	}
	if(onl.chn < 0 || 127 < onl.chn)
	{
		return false;
	}

	return true;
}

bool InttMap::IsValidOrWildcard(
	InttMap::Online_s const& onl)
{
	if(onl.lyr != InttMap::Wildcard && (onl.lyr < 0 || 3 < onl.lyr))
	{
		return false;
	}
	if(onl.ldr != InttMap::Wildcard && (onl.ldr < 0 || (onl.lyr < 2 ? 11 : 15) < onl.ldr))
	{
		return false;
	}
	if(onl.arm != InttMap::Wildcard && (onl.arm < 0 || 1 < onl.arm))
	{
		return false;
	}
	if(onl.chp != InttMap::Wildcard && (onl.chp < 0 || 25 < onl.chp))
	{
		return false;
	}
	if(onl.chn != InttMap::Wildcard && (onl.chn < 0 || 127 < onl.chn))
	{
		return false;
	}

	return true;
}

const struct InttMap::Online_s InttMap::OnlineBegin {
	.lyr = 0,
	.ldr = 0,
	.arm = 0,
	.chp = 0,
	.chn = 0,
};

const struct InttMap::Online_s InttMap::OnlineEnd {
	.lyr = 4,
	.ldr = 16,
	.arm = 4,
	.chp = 8,
	.chn = 256,
};

const struct InttMap::Online_s InttMap::OnlineRBegin {
	.lyr = 3,
	.ldr = 15,
	.arm = 1,
	.chp = 25,
	.chn = 127,
};

const struct InttMap::Online_s InttMap::OnlineREnd {
	.lyr = 4,
	.ldr = 16,
	.arm = 4,
	.chp = 8,
	.chn = 256,
};

InttMap::RawData_s& InttMap::RawData_s::operator++()
{
	if(!InttMap::IsValidOrWildcard(*this))
	{
		*this = InttMap::RawDataEnd;
		return *this;
	}

	if(chn != InttMap::Wildcard)
	{
		if(++chn < 128)
		{
			return *this;
		}
		chn = 0;
	}
	if(chp != InttMap::Wildcard)
	{
		if(++chp < 26)
		{
			return *this;
		}
		chp = 0;
	}
	if(fee != InttMap::Wildcard)
	{
		if(++fee < 14)
		{
			return *this;
		}
		fee = 0;
	}
	if(pid != InttMap::Wildcard)
	{
		if(++pid < 3009)
		{
			return *this;
		}
		pid = 3001;
	}

	*this = InttMap::RawDataEnd;
	return *this;
}

InttMap::RawData_s& InttMap::RawData_s::operator--()
{
	if(!InttMap::IsValidOrWildcard(*this))
	{
		*this = InttMap::RawDataREnd;
		return *this;
	}

	if(chn != InttMap::Wildcard)
	{
		if(--chn >= 0)
		{
			return *this;
		}
		chn = 127;
	}
	if(chp != InttMap::Wildcard)
	{
		if(--chp >= 0)
		{
			return *this;
		}
		chp = 25;
	}
	if(fee != InttMap::Wildcard)
	{
		if(--fee >= 0)
		{
			return *this;
		}
		fee = 13;
	}
	if(pid != InttMap::Wildcard)
	{
		if(--pid >= 3001)
		{
			return *this;
		}
		pid = 3008;
	}

	*this = InttMap::RawDataREnd;
	return *this;
}

bool InttMap::IsValid(
	InttMap::RawData_s const& raw)
{
	if(raw.pid < 3001 || 3008 < raw.pid)
	{
		return false;
	}
	if(raw.fee < 0 || 13 < raw.fee)
	{
		return false;
	}
	if(raw.chp < 0 || 25 < raw.chp)
	{
		return false;
	}
	if(raw.chn < 0 || 127 < raw.chn)
	{
		return false;
	}

	return true;
}

bool InttMap::IsValidOrWildcard(
	InttMap::RawData_s const& raw)
{
	if(raw.pid != InttMap::Wildcard && (raw.pid < 3001 || 3008 < raw.pid))
	{
		return false;
	}
	if(raw.fee != InttMap::Wildcard && (raw.fee < 0 || 13 < raw.fee))
	{
		return false;
	}
	if(raw.chp != InttMap::Wildcard && (raw.chp < 0 || 25 < raw.chp))
	{
		return false;
	}
	if(raw.chn != InttMap::Wildcard && (raw.chn < 0 || 127 < raw.chn))
	{
		return false;
	}

	return true;
}

const struct InttMap::RawData_s InttMap::RawDataBegin {
	.pid = 3001,
	.fee = 0,
	.chp = 0,
	.chn = 0,
};

const struct InttMap::RawData_s InttMap::RawDataEnd {
	.pid = 3009,
	.fee = 14,
	.chp = 26,
	.chn = 256,
};

const struct InttMap::RawData_s InttMap::RawDataRBegin {
	.pid = 3008,
	.fee = 13,
	.chp = 25,
	.chn = 127,
};

const struct InttMap::RawData_s InttMap::RawDataREnd {
	.pid = 3009,
	.fee = 14,
	.chp = 26,
	.chn = 256,
};

InttMap::Offline_s& InttMap::Offline_s::operator++()
{
	if(!InttMap::IsValidOrWildcard(*this))
	{
		*this = InttMap::OfflineEnd;
		return *this;
	}

	if(strip_phi != InttMap::Wildcard)
	{
		if(++strip_phi < 256)
		{
			return *this;
		}
		strip_phi = 0;
	}
	if(strip_z != InttMap::Wildcard)
	{
		if(++strip_z < (ladder_z % 2 ? 5 : 8))
		{
			return *this;
		}
		strip_z = 0;
	}
	if(ladder_z != InttMap::Wildcard)
	{
		if(++ladder_z < 4)
		{
			return *this;
		}
		ladder_z = 0;
	}
	if(ladder_phi != InttMap::Wildcard)
	{
		if(++ladder_phi < (layer < 5 ? 12 : 16))
		{
			return *this;
		}
		ladder_phi = 0;
	}
	if(layer != InttMap::Wildcard)
	{
		if(++layer < 7)
		{
			return *this;
		}
		layer = 3;
	}

	*this = InttMap::OfflineEnd;
	return *this;
}

InttMap::Offline_s& InttMap::Offline_s::operator--()
{
	if(!InttMap::IsValidOrWildcard(*this))
	{
		*this = InttMap::OfflineREnd;
		return *this;
	}

	if(strip_phi != InttMap::Wildcard)
	{
		if(--strip_phi >= 0)
		{
			return *this;
		}
		strip_phi = 255;
	}
	if(strip_z != InttMap::Wildcard)
	{
		if(--strip_z >= 0)
		{
			return *this;
		}
		strip_z = ladder_z % 2 ? 7 : 4;
	}
	if(ladder_z != InttMap::Wildcard)
	{
		if(--ladder_z >= 0)
		{
			return *this;
		}
		ladder_z = 3;
	}
	if(ladder_phi != InttMap::Wildcard)
	{
		if(--ladder_phi >= 0)
		{
			return *this;
		}
		ladder_phi = layer < 6 ? 11 : 15;
	}
	if(layer != InttMap::Wildcard)
	{
		if(--layer >= 3)
		{
			return *this;
		}
		layer = 6;
	}

	*this = InttMap::OfflineREnd;
	return *this;
}

bool InttMap::IsValid(
	InttMap::Offline_s const& ofl)
{
	if(ofl.layer < 3 || 6 < ofl.layer)
	{
		return false;
	}
	if(ofl.ladder_phi < 0 || (ofl.layer < 5 ? 11 : 15) < ofl.ladder_phi)
	{
		return false;
	}
	if(ofl.ladder_z < 0 || 3 < ofl.ladder_z)
	{
		return false;
	}
	if(ofl.strip_z < 0 || (ofl.ladder_z % 2 ? 4 : 7) < ofl.strip_z)
	{
		return false;
	}
	if(ofl.strip_phi < 0 || 255 < ofl.strip_phi)
	{
		return false;
	}

	return true;
}

bool InttMap::IsValidOrWildcard(
	InttMap::Offline_s const& ofl)
{
	if(ofl.layer != InttMap::Wildcard && (ofl.layer < 3 || 6 < ofl.layer))
	{
		return false;
	}
	if(ofl.ladder_phi != InttMap::Wildcard && (ofl.ladder_phi < 0 || (ofl.layer < 5 ? 11 : 15) < ofl.ladder_phi))
	{
		return false;
	}
	if(ofl.ladder_z != InttMap::Wildcard && (ofl.ladder_z < 0 || 3 < ofl.ladder_z))
	{
		return false;
	}
	if(ofl.strip_z != InttMap::Wildcard && (ofl.strip_z < 0 || (ofl.ladder_z % 2 ? 4 : 7) < ofl.strip_z))
	{
		return false;
	}
	if(ofl.strip_phi != InttMap::Wildcard && (ofl.strip_phi < 0 || 255 < ofl.strip_phi))
	{
		return false;
	}

	return true;
}

const struct InttMap::Offline_s InttMap::OfflineBegin {
	.layer = 3,
	.ladder_phi = 0,
	.ladder_z = 0,
	.strip_z = 0,
	.strip_phi = 0,
};

const struct InttMap::Offline_s InttMap::OfflineEnd {
	.layer = 7,
	.ladder_phi = 16,
	.ladder_z = 4,
	.strip_z = 8,
	.strip_phi = 256,
};

const struct InttMap::Offline_s InttMap::OfflineRBegin {
	.layer = 6,
	.ladder_phi = 15,
	.ladder_z = 3,
	.strip_z = 4,
	.strip_phi = 255,
};

const struct InttMap::Offline_s InttMap::OfflineREnd {
	.layer = 7,
	.ladder_phi = 16,
	.ladder_z = 4,
	.strip_z = 8,
	.strip_phi = 256,
};

bool InttMap::OnlineComparator::operator()(
    struct InttMap::Online_s const& lhs,
    struct InttMap::Online_s const& rhs) const
{
  if (lhs.lyr != rhs.lyr)
  {
    return lhs.lyr < rhs.lyr;
  }
  if (lhs.ldr != rhs.ldr)
  {
    return lhs.ldr < rhs.ldr;
  }
  if (lhs.arm != rhs.arm)
  {
    return lhs.arm < rhs.arm;
  }
  if (lhs.chp != rhs.chp)
  {
    return lhs.chp < rhs.chp;
  }
  if (lhs.chn != rhs.chn)
  {
    return lhs.chn < rhs.chn;
  }

  return false;
}

bool InttMap::OnlineWildcardComparator::operator()(
    struct InttMap::Online_s const& lhs,
    struct InttMap::Online_s const& rhs) const
{
  if (lhs.lyr == Wildcard || rhs.lyr == Wildcard)
  {
    return false;
  }
  if (lhs.lyr != rhs.lyr)
  {
    return lhs.lyr < rhs.lyr;
  }

  if (lhs.ldr == Wildcard || rhs.ldr == Wildcard)
  {
    return false;
  }
  if (lhs.ldr != rhs.ldr)
  {
    return lhs.ldr < rhs.ldr;
  }

  if (lhs.arm == Wildcard || rhs.arm == Wildcard)
  {
    return false;
  }
  if (lhs.arm != rhs.arm)
  {
    return lhs.arm < rhs.arm;
  }

  if (lhs.chp == Wildcard || rhs.chp == Wildcard)
  {
    return false;
  }
  if (lhs.chp != rhs.chp)
  {
    return lhs.chp < rhs.chp;
  }

  if (lhs.chn == Wildcard || rhs.chn == Wildcard)
  {
    return false;
  }
  if (lhs.chn != rhs.chn)
  {
    return lhs.chn < rhs.chn;
  }

  return false;
}

bool InttMap::RawDataComparator::operator()(
    struct InttMap::RawData_s const& lhs,
    struct InttMap::RawData_s const& rhs) const
{
  if (lhs.pid != rhs.pid)
  {
    return lhs.pid < rhs.pid;
  }
  if (lhs.fee != rhs.fee)
  {
    return lhs.fee < rhs.fee;
  }
  if (lhs.chp != rhs.chp)
  {
    return lhs.chp < rhs.chp;
  }
  if (lhs.chn != rhs.chn)
  {
    return lhs.chn < rhs.chn;
  }

  return false;
}

bool InttMap::RawDataWildcardComparator::operator()(
    struct InttMap::RawData_s const& lhs,
    struct InttMap::RawData_s const& rhs) const
{
  if (lhs.pid == Wildcard || rhs.pid == Wildcard)
  {
    return false;
  }
  if (lhs.pid != rhs.pid)
  {
    return lhs.pid < rhs.pid;
  }

  if (lhs.fee == Wildcard || rhs.fee == Wildcard)
  {
    return false;
  }
  if (lhs.fee != rhs.fee)
  {
    return lhs.fee < rhs.fee;
  }

  if (lhs.chp == Wildcard || rhs.chp == Wildcard)
  {
    return false;
  }
  if (lhs.chp != rhs.chp)
  {
    return lhs.chp < rhs.chp;
  }

  if (lhs.chn == Wildcard || rhs.chn == Wildcard)
  {
    return false;
  }
  if (lhs.chn != rhs.chn)
  {
    return lhs.chn < rhs.chn;
  }

  return false;
}

bool InttMap::OfflineComparator::operator()(
    struct InttMap::Offline_s const& lhs,
    struct InttMap::Offline_s const& rhs) const
{
  if (lhs.layer != rhs.layer)
  {
    return lhs.layer < rhs.layer;
  }
  if (lhs.ladder_phi != rhs.ladder_phi)
  {
    return lhs.ladder_phi < rhs.ladder_phi;
  }
  if (lhs.ladder_z != rhs.ladder_z)
  {
    return lhs.ladder_z < rhs.ladder_z;
  }
  if (lhs.strip_z != rhs.strip_z)
  {
    return lhs.strip_z < rhs.strip_z;
  }
  if (lhs.strip_phi != rhs.strip_phi)
  {
    return lhs.strip_phi < rhs.strip_phi;
  }

  return false;
}

bool InttMap::OfflineWildcardComparator::operator()(
    struct InttMap::Offline_s const& lhs,
    struct InttMap::Offline_s const& rhs) const
{
  if (lhs.layer == Wildcard || rhs.layer == Wildcard)
  {
    return false;
  }
  if (lhs.layer != rhs.layer)
  {
    return lhs.layer < rhs.layer;
  }

  if (lhs.ladder_phi == Wildcard || rhs.ladder_phi == Wildcard)
  {
    return false;
  }
  if (lhs.ladder_phi != rhs.ladder_phi)
  {
    return lhs.ladder_phi < rhs.ladder_phi;
  }

  if (lhs.ladder_z == Wildcard || rhs.ladder_z == Wildcard)
  {
    return false;
  }
  if (lhs.ladder_z != rhs.ladder_z)
  {
    return lhs.ladder_z < rhs.ladder_z;
  }

  if (lhs.strip_z == Wildcard || rhs.strip_z == Wildcard)
  {
    return false;
  }
  if (lhs.strip_z != rhs.strip_z)
  {
    return lhs.strip_z < rhs.strip_z;
  }

  if (lhs.strip_phi == Wildcard || rhs.strip_phi == Wildcard)
  {
    return false;
  }
  if (lhs.strip_phi != rhs.strip_phi)
  {
    return lhs.strip_phi < rhs.strip_phi;
  }

  return false;
}

bool operator<(
  struct InttMap::Online_s const& lhs,
  struct InttMap::Online_s const& rhs)
{
  if(lhs.lyr != rhs.lyr)
  {
	return lhs.lyr < rhs.lyr;
  }
  if(lhs.ldr != rhs.ldr)
  {
    return lhs.ldr < rhs.ldr;
  }
  if(lhs.arm != rhs.arm)
  {
    return lhs.arm < rhs.arm;
  }
  if(lhs.chp != rhs.chp)
  {
    return lhs.chp < rhs.chp;
  }
  if(lhs.chn != rhs.chn)
  {
    return lhs.chn < rhs.chn;
  }

  return false;
}

bool operator==(
  struct InttMap::Online_s const& lhs,
  struct InttMap::Online_s const& rhs)
{
  if(lhs.lyr != rhs.lyr)
  {
    return false;
  }
  if(lhs.ldr != rhs.ldr)
  {
    return false;
  }
  if(lhs.arm != rhs.arm)
  {
    return false;
  }
  if(lhs.chp != rhs.chp)
  {
    return false;
  }
  if(lhs.chn != rhs.chn)
  {
    return false;
  }

  return true;
}

bool operator!=(
  struct InttMap::Online_s const& lhs,
  struct InttMap::Online_s const& rhs)
{
  return !(lhs == rhs);
}

bool operator>(
  struct InttMap::Online_s const& lhs,
  struct InttMap::Online_s const& rhs)
{
  return rhs < lhs;
}

bool operator<=(
  struct InttMap::Online_s const& lhs,
  struct InttMap::Online_s const& rhs)
{
  return !(lhs > rhs);
}

bool operator>=(
  struct InttMap::Online_s const& lhs,
  struct InttMap::Online_s const& rhs)
{
  return !(lhs < rhs);
}

bool operator<(
  struct InttMap::RawData_s const& lhs,
  struct InttMap::RawData_s const& rhs)
{
  if(lhs.pid != rhs.pid)
  {
    return lhs.pid < rhs.pid;
  }
  if(lhs.fee != rhs.fee)
  {
    return lhs.fee < rhs.fee;
  }
  if(lhs.chp != rhs.chp)
  {
    return lhs.chp < rhs.chp;
  }
  if(lhs.chn != rhs.chn)
  {
    return lhs.chn < rhs.chn;
  }

  return false;
}

bool operator==(
  struct InttMap::RawData_s const& lhs,
  struct InttMap::RawData_s const& rhs)
{
  if(lhs.pid != rhs.pid)
  {
    return false;
  }
  if(lhs.fee != rhs.fee)
  {
    return false;
  }
  if(lhs.chp != rhs.chp)
  {
    return false;
  }
  if(lhs.chn != rhs.chn)
  {
    return false;
  }

  return true;
}

bool operator!=(
  struct InttMap::RawData_s const& lhs,
  struct InttMap::RawData_s const& rhs)
{
  return !(lhs == rhs);
}

bool operator>(
  struct InttMap::RawData_s const& lhs,
  struct InttMap::RawData_s const& rhs)
{
  return rhs < lhs;
}

bool operator<=(
  struct InttMap::RawData_s const& lhs,
  struct InttMap::RawData_s const& rhs)
{
  return !(lhs > rhs);
}

bool operator>=(
  struct InttMap::RawData_s const& lhs,
  struct InttMap::RawData_s const& rhs)
{
  return !(lhs < rhs);
}

bool operator<(
  struct InttMap::Offline_s const& lhs,
  struct InttMap::Offline_s const& rhs)
{
  if(lhs.layer != rhs.layer)
  {
    return lhs.layer < rhs.layer;
  }
  if(lhs.ladder_phi != rhs.ladder_phi)
  {
    return lhs.ladder_phi < rhs.ladder_phi;
  }
  if(lhs.ladder_z != rhs.ladder_z)
  {
    return lhs.ladder_z < rhs.ladder_z;
  }
  if(lhs.strip_z != rhs.strip_z)
  {
    return lhs.strip_z < rhs.strip_z;
  }
  if(lhs.strip_phi != rhs.strip_phi)
  {
    return lhs.strip_phi < rhs.strip_phi;
  }

  return false;
}

bool operator==(
  struct InttMap::Offline_s const& lhs,
  struct InttMap::Offline_s const& rhs)
{
  if(lhs.layer != rhs.layer)
  {
    return false;
  }
  if(lhs.ladder_phi != rhs.ladder_phi)
  {
    return false;
  }
  if(lhs.ladder_z != rhs.ladder_z)
  {
    return false;
  }
  if(lhs.strip_z != rhs.strip_z)
  {
    return false;
  }
  if(lhs.strip_phi != rhs.strip_phi)
  {
    return false;
  }

  return true;
}

bool operator!=(
  struct InttMap::Offline_s const& lhs,
  struct InttMap::Offline_s const& rhs)
{
  return !(lhs == rhs);
}

bool operator>(
  struct InttMap::Offline_s const& lhs,
  struct InttMap::Offline_s const& rhs)
{
  return rhs < lhs;
}

bool operator<=(
  struct InttMap::Offline_s const& lhs,
  struct InttMap::Offline_s const& rhs)
{
  return !(lhs > rhs);
}

bool operator>=(
  struct InttMap::Offline_s const& lhs,
  struct InttMap::Offline_s const& rhs)
{
  return !(lhs < rhs);
}

std::ostream& operator<<(
  std::ostream& os,
  struct InttMap::Online_s const& onl)
{
  os << "InttMap::Online_s\n"
	 << boost::str(boost::format("\t%-12s%4d\n") % "lyr:" % onl.lyr)
	 << boost::str(boost::format("\t%-12s%4d\n") % "ldr:" % onl.ldr)
	 << boost::str(boost::format("\t%-12s%4d\n") % "arm:" % onl.arm)
	 << boost::str(boost::format("\t%-12s%4d\n") % "chp:" % onl.chp)
	 << boost::str(boost::format("\t%-12s%4d\n") % "chn:" % onl.chn);
  return os;
}

std::ostream& operator<<(
  std::ostream& os,
  struct InttMap::RawData_s const& raw)
{
  os << "InttMap::RawData_s\n"
	 << boost::str(boost::format("\t%-12s%4d\n") % "pid:" % raw.pid)
	 << boost::str(boost::format("\t%-12s%4d\n") % "fee:" % raw.fee)
	 << boost::str(boost::format("\t%-12s%4d\n") % "chp:" % raw.chp)
	 << boost::str(boost::format("\t%-12s%4d\n") % "chn:" % raw.chn);
  return os;
}

std::ostream& operator<<(
  std::ostream& os,
  struct InttMap::Offline_s const& ofl)
{
  os << "InttMap::Offline_s\n"
	 << boost::str(boost::format("\t%-12s%4d\n") % "layer:" % ofl.layer)
	 << boost::str(boost::format("\t%-12s%4d\n") % "ladder_phi:" % ofl.ladder_phi)
	 << boost::str(boost::format("\t%-12s%4d\n") % "ladder_z:" % ofl.ladder_z)
	 << boost::str(boost::format("\t%-12s%4d\n") % "strip_z:" % ofl.strip_z)
	 << boost::str(boost::format("\t%-12s%4d\n") % "strip_phi:" % ofl.strip_phi);
  return os;
}
