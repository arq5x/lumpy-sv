#ifndef _XAMTAG_H
#define _XAMTAG_H

#include <string>

/**
 * @brief class to represent a Xam TAG:TYPE:VALUE entry
 */
template<class TAG_TYPE>
class XamTag {
    public:
        explicit XamTag(const std::string& name, const TAG_TYPE& value, const bool missing = false) :
            m_name(name),
            m_value(value),
            m_missing(missing)
        {}

        explicit XamTag(const std::string& name, TAG_TYPE&& value, const bool missing = false) :
            m_name(name),
            m_value(std::move(value)),
            m_missing(missing)
        {}

        XamTag(const XamTag& other) = default;
        XamTag(XamTag&& other) = default;
        XamTag& operator=(const XamTag& other) = default;
        XamTag& operator=(XamTag&& other) = default;
        ~XamTag() = default;

        std::string name() const { return m_name; }
        TAG_TYPE value() const { return m_value; }
        bool missing() const { return m_missing; }

    private:
        std::string m_name;
        TAG_TYPE m_value;
        bool m_missing;
};


#endif 